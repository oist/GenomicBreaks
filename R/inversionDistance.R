#' Inversion distance
#'
#' Implementation of Hannenhalli and Pevzner algorithm
#'
#' @param p A permutation vector
#'
#' @returns The minimal number of inversions required to sort the given permutation vector
#'
#' @export

# Defining all required functions for the inversion distance

#install.packages("igraph")
#library(igraph)





##### Extended Permutation #####

extendedPermutation <- function(p){
  p_extended <- sapply(p, function(i) {
    if (i>0) { c( 2L * i,        2L * i + 1L) }
    else     { c(-2L * i + 1L , -2L * i     ) }
  })
  c(1L, p_extended, 2L * length(p) + 2L)
}





##### Breakpoint Graph #####

breakpoint_graph <- function(p_extended) {
  len_pe <- length(p_extended)
  len_p <- len_pe / 2 - 1

  # initialize edge lists
  edge_list <- c()
  edge_colors <- c()
  edge_unoriented <- c()

  for (i in 0:len_p) {
    a <- p_extended[2 * i + 1]
    b <- p_extended[2 * i + 2]

    # black edge if not consecutive in genome permutation
    if (abs(a - b) != 1) {
      edge_list <- c(edge_list, a, b)
      edge_colors <- c(edge_colors, "black")
      edge_unoriented <- c(edge_unoriented, 1)
    }

    # gray edge if not adjacent in the list
    i_init <- match(2 * i + 1, p_extended)
    i_end  <- match(2 * i + 2, p_extended)

    if (abs(i_end - i_init) != 1) {
      edge_list <- c(edge_list, 2 * i + 1, 2 * i + 2)
      edge_colors <- c(edge_colors, "gray")
      edge_unoriented <- c(edge_unoriented, (i_init + i_end) %% 2)
    }
  }

  # build graph and assign attributes
  g <- make_empty_graph(n = len_pe, directed = FALSE)
  g <- add_edges(g, edge_list)
  E(g)$color <- edge_colors
  E(g)$unoriented <- edge_unoriented

  return(g)
}



# Check if two edges interleaves

is_interleaving <- function(p_extended, graph_1, graph_2, edge_1, edge_2) {
  # get positions of the two ends of each edge in the extended permutation
  e1_nodes <- ends(graph_1, edge_1)
  e2_nodes <- ends(graph_2, edge_2)

  e1_pos <- sort(match(e1_nodes, p_extended))
  e2_pos <- sort(match(e2_nodes, p_extended))

  # check for interleaving: overlap without containment
  return((e1_pos[1] < e2_pos[1] && e2_pos[1] < e1_pos[2] && e1_pos[2] < e2_pos[2]) ||
           (e2_pos[1] < e1_pos[1] && e1_pos[1] < e2_pos[2] && e2_pos[2] < e1_pos[2]))
}



# Count the number of breakpoints (including artificial extremities)

bp_count <- function(p_extended){
  sum(abs(diff(p_extended)) != 1 )}



# Count the number of cycles

cycle_count <- function(g){
  sum(components(g)$csize > 1)
}





##### Components Graph #####

# Connected components graph given the breakpoint graph

components_graph <- function(g, p_extended) {
  comps <- components(g)
  comp_ids <- which(comps$csize > 1)

  # map components to vertices
  component_vertices <- lapply(comp_ids, function(i) which(comps$membership == i))
  component_subgraphs <- lapply(component_vertices, function(vids) induced_subgraph(g, vids))

  unoriented_flags <- sapply(component_subgraphs, function(cycle) all(E(cycle)$unoriented == 1))
  colors <- ifelse(unoriented_flags, "lightblue", "lightsalmon")
  min_extents <- sapply(component_vertices, function(verts) min(match(verts, p_extended)))
  max_extents <- sapply(component_vertices, function(verts) max(match(verts, p_extended)))

  # build h with attributes
  h <- make_empty_graph(n = length(comp_ids), directed = FALSE)
  V(h)$cycle_index <- comp_ids
  V(h)$unoriented <- unoriented_flags
  V(h)$color <- colors
  V(h)$min_extent <- min_extents
  V(h)$max_extent <- max_extents

  if (length(h) > 1) {
    for (i in 1:(length(h) - 1)) {
      vertices_i <- component_vertices[[i]]
      cycle_i <- component_subgraphs[[i]]
      V(cycle_i)$name <- vertices_i  # set names for is_interleaving

      gray_i <- E(cycle_i)[E(cycle_i)$color == "gray"]

      for (j in (i + 1):length(h)) {
        vertices_j <- component_vertices[[j]]
        cycle_j <- component_subgraphs[[j]]
        V(cycle_j)$name <- vertices_j  # same

        gray_j <- E(cycle_j)[E(cycle_j)$color == "gray"]

        found <- FALSE
        for (e1 in gray_i) {
          for (e2 in gray_j) {
            if (is_interleaving(p_extended, cycle_i, cycle_j, e1, e2)) {
              h <- add_edges(h, c(i, j))
              found <- TRUE
              break
            }
          }
          if (found) break
        }
      }
    }
  }


  return(h)
}





##### Identifying hurdles #####

hurdles_count <- function(g, p_extended) {
  h <- components_graph(g, p_extended)

  comps_h <- components(h)
  comps_g <- components(g)
  n_comps_h <- max(comps_h$membership)

  # Identify all unoriented components in h
  is_unoriented <- sapply(1:n_comps_h, function(i) {
    all(V(h)$unoriented[comps_h$membership == i])
  })
  unoriented_ids <- which(is_unoriented)

  n_unoriented <- length(unoriented_ids)
  if (n_unoriented == 0) {
    return(data.frame(
      membership = integer(0), min_extent = numeric(0), max_extent = numeric(0),
      contained = integer(0), contains = integer(0),
      hurdle = integer(0), superhurdle = integer(0)
    ))
  }

  # Preallocate info table
  info <- data.frame(
    membership = unoriented_ids,
    min_extent = numeric(n_unoriented),
    max_extent = numeric(n_unoriented),
    contained = integer(n_unoriented),
    contains = integer(n_unoriented),
    hurdle = integer(n_unoriented),
    superhurdle = integer(n_unoriented)
  )

  # Fill min and max extents
  for (i in seq_along(unoriented_ids)) {
    verts <- which(comps_h$membership == unoriented_ids[i])
    info$min_extent[i] <- min(V(h)$min_extent[verts])
    info$max_extent[i] <- max(V(h)$max_extent[verts])
  }

  if (n_unoriented == 1) {
    info$hurdle[1] <- 1
    return(info)
  }

  # Determine containment relationships
  contains_matrix <- outer(info$min_extent, info$min_extent, "<") &
    outer(info$max_extent, info$max_extent, ">")
  diag(contains_matrix) <- FALSE

  info$contains <- rowSums(contains_matrix)
  info$contained <- colSums(contains_matrix)

  # Determine hurdles
  for (i in seq_len(n_unoriented)) {
    if ((info$contained[i] == 0 || info$contained[i] == 1) && info$contains[i] == 0) {
      info$hurdle[i] <- 1
    }
    if (info$contained[i] == 0 && info$contains[i] == 1) {
      info$hurdle[i] <- 2  # possible superhurdle
    }
  }

  if (any(info$hurdle == 2)) {
    # Check for actual superhurdle status
    info_smallest <- which(info$hurdle == 1)
    min_hurdles <- min(info$min_extent[info_smallest])
    max_hurdles <- max(info$max_extent[info_smallest])

    possible_greatest <- which(info$hurdle == 2)
    for (i in possible_greatest) {
      if (info$min_extent[i] > min_hurdles || info$max_extent[i] < max_hurdles) {
        info$hurdle[i] <- 0
      }
    }

    # If still a potential superhurdle, validate separation
    super_candidate <- which(info$hurdle == 2)
    if (length(super_candidate) == 1) {
      sc_idx <- super_candidate[1]
      info$hurdle[sc_idx] <- 1  # start assuming it's valid

      comp_id <- info$membership[sc_idx]
      comp_vertices <- which(comps_h$membership == comp_id)
      cycle_ids <- V(h)$cycle_index[comp_vertices]

      if (length(info_smallest) > 1) {
        edge_colors <- E(g)$color == "gray"
        gray_edges <- which(edge_colors)
        edge_ends <- ends(g, gray_edges, names = FALSE)
        positions <- match(edge_ends, p_extended)

        for (i in cycle_ids) {
          edge_indices <- which(comps_g$membership[gray_edges] == i)
          if (length(edge_indices) == 0) next

          for (j in 1:(length(info_smallest)-1)) {
            for (k in (j+1):length(info_smallest)) {
              min_j <- info$min_extent[info_smallest[j]]
              max_j <- info$max_extent[info_smallest[j]]
              min_k <- info$min_extent[info_smallest[k]]
              max_k <- info$max_extent[info_smallest[k]]

              for (e_idx in edge_indices) {
                edge_pos <- sort(positions[e_idx, ])

                condition_1 <- (edge_pos[1] < min_j && edge_pos[2] > max_j &&
                                  (edge_pos[1] > max_k || edge_pos[2] < min_k))
                condition_2 <- (edge_pos[1] < min_k && edge_pos[2] > max_k &&
                                  (edge_pos[1] > max_j || edge_pos[2] < min_j))

                if (condition_1 || condition_2) {
                  info$hurdle[sc_idx] <- 0
                  break
                }
              }

              if (info$hurdle[sc_idx] == 0) break
            }
            if (info$hurdle[sc_idx] == 0) break
          }
          if (info$hurdle[sc_idx] == 0) break
        }
      }
    }
  }

  return(info)
}





##### Identifying superhurdles #####

superhurdles_count <- function(info, g, query_sequence_unsig){

  #if the number of hurdles is even, it is not necessary to calculate superhurdles because of the definition of fortress

  if (sum(info$hurdle)%%2 == 0) {return(info)}

  #if all unoriented component are hurdles, there is no hurdle protecting an unoriented nonhurdle

  if (sum(info$hurdle) == length(info$hurdle)) {return(info)}

  #for all hurdles: remove one by one and check if some non hurdle turned into a hurdle; if it is the case, then the removed hurdle is a superhurdle

  for (ignored in info$membership[info$hurdle == 1]){

    # remove "ignored" to calculate hurdles

    info2 <- subset(info, membership != ignored)

    line <- length(info2$membership)

    if (line == 1){
      info["superhurdle"][info["membership"] == ignored] <- 1
      return(info)
    }

    for (i in 1:line){
      info2[i, c("contained", "contains", "hurdle", "superhurdle")] <- c(0, 0, 0, 0)
    }

    for (i in 1:(line-1)) {
      for (j in 2:line) {

        min_i <- info2[i, "min_extent"]
        max_i <- info2[i, "max_extent"]

        min_j <- info2[j, "min_extent"]
        max_j <- info2[j, "max_extent"]

        if ((min_i < min_j) & (min_j < max_j) & (max_j < max_i)) {
          info2[i, "contains"] <- 1
          info2[j, "contained"] <- 1
        }

        if ((min_j < min_i) & (min_i < max_i) & (max_i< max_j)) {
          info2[j, "contains"] <- 1
          info2[i, "contained"] <- 1
        }
      }
    }

    #hurdles never are in "the middle". They need to contain and not being contained,
    # or being contained and not contain

    for (i in 1:line) {
      if ((info2[i, "contained"] == 0 || info2[i, "contained"] == 1) && info2[i, "contains"] == 0){
        info2[i, "hurdle"] <- 1
      }
      if (info2[i, "contained"] == 0 && info2[i, "contains"] == 1){
        info2[i, "hurdle"] <- 2 #we need to analyze this case
      }
    }

    if (any(info2$hurdle == 2)) {

      #first, check if it covers all other hurdles, otherwise it is not a (greatest) hurdle

      possible_greatest <- which(info2$hurdle == 2)

      min_hurdles <-  min(info2[info2$hurdle == 1, ]["min_extent"])
      max_hurdles <-  max(info2[info2$hurdle == 1, ]["max_extent"])

      for (i in possible_greatest){

        if (info2[i, "min_extent"] > min_hurdles || info2[i, "max_extent"] < max_hurdles){
          info2[i, "hurdle"] <- 0
        }
      }

      #if there is still someone where hurdle == 2 (it is unique), check the steps below

      if (any(info2$hurdle == 2)) {

        #attribute cycle_index from vertices in h allows to access the correspondent cycle in g

        greatest <- which(info2$hurdle == 2)
        greatest <- info2[greatest, "membership"]

        comp_vertices <- which(comps$membership == greatest)

        idx = V(h)$cycle_index[comp_vertices] #get the memberships of the cycles which compose that component

        info_smallest <- which(info2$hurdle == 1)

        #we start assuming it is a hurdle; if we find it dividing two hurdles, we put a zero.

        info2[greatest, "hurdle"] <- 1

        #if info_smallest has only 1 hurdle, it means the great can't separate two hurdles

        if (length(info_smallest) == 1){
          return(info)
        }

        greatest_check <- FALSE

        for (i in idx) {

          if(!(greatest_check)){

            #get edges of that component
            edges_to_analyze <- E(g)[components(g)$membership == i & E(g)$color=="gray"]

            #get all pairs of another hurdles

            for (j in info_smallest[1:(length(info_smallest)-1)]) {
              for (k in info_smallest[2:(length(info_smallest))]) {

                min_j <- info2[j, "min_extent"]  # extent(U')[1]
                max_j <- info2[j, "max_extent"]  # extent(U')[2]

                min_k <- info2[k, "min_extent"]  # extent(U")[1]
                max_k <- info2[k, "max_extent"]  # extent(U")[2]

                for (edge in edges_to_analyze){

                  # gray edge from U
                  edge_extent <- sort(match(c(ends(g, edge)[1], ends(g, edge)[2]), query_sequence_unsig))

                  # U separates U' and U'' if there is a gray edge in U containing extent(U'),
                  # but without intersection with extent(U") (and vice versa)

                  condition_1 <- ((edge_extent[1] < min_j & edge_extent[2] > max_j) &
                                    (edge_extent[1] > max_k || edge_extent[2] < min_k))

                  condition_2 <- ((edge_extent[1] < min_k & edge_extent[2] > max_k) &
                                    (edge_extent[1] > max_j || edge_extent[2] < min_j))

                  if (condition_1 || condition_2){
                    info2[greatest, "hurdle"] <- 0
                    greatest_check <- TRUE
                  }
                }
              }
            }
          }
        }
      }
    }

    #compare info and info2 to see if some nonhurdle turned into a hurdle

    if (any(((subset(info, membership != ignored)$hurdle) - (info2$hurdle)) == -1)){

      #ignored is a superhurdle

      info["superhurdle"][info["membership"] == ignored] <- 1
    }
  }

  return(info)
}

#analyze if the permutation is a fortress
is_fortress <- function(superhurdles){

  if(sum(superhurdles$hurdles)%%2 == 0){
    return(0)
  }
  if(sum(superhurdles$superhurdle) == sum(superhurdles$hurdle)){
    return(1)
  }
  return(0)

  #a permutation is a fortress if the number of hurdles is odd and all hurdles are superhurdles
  #return 1 if it is fortress
  #return 0 otherwise

}



##### Inversion Distance #####

inversionDistance <- function(p){

  p_extended <- extendedPermutation(p)

  bp_graph <- breakpoint_graph(p_extended)
  hurdles <- hurdles_count(bp_graph, p_extended)
  superhurdles <- superhurdles_count(hurdles, bp_graph, p_extended)

  minimal <- bp_count(p_extended) - cycle_count(bp_graph) + sum(hurdles$hurdle) + is_fortress(superhurdles)

  return(minimal)
}
