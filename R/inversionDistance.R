#' Inversion Distance
#'
#' Computes the minimal number of inversions required to sort a signed permutation
#' using the Hannenhalli and Pevzner algorithm.
#'
#' This function uses several internal helper functions (e.g., `extendedPermutation`,
#' `breakpoint_graph`, `hurdles_count`, `superhurdles_count`, and others) to
#' compute properties of the breakpoint graph and identify cycles, hurdles and superhurdles.
#' It also depends on `permutationVector()`, which is defined in another .R file.
#'
#' @param gb_object A GBreaks object containing genome data from which a permutation vector can be extracted
#' using the \code{permutationVector()} function.
#'
#' @return An integer: the minimal number of inversions needed to sort the permutation.
#'
#' @seealso \code{\link{permutationVector}} for generating the permutation vector.
#'
#' @import igraph
#'
#' @author Bruna Fistarol
#'
#' @export





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

  # Initialize edge lists
  edge_list <- c()
  edge_colors <- c()
  edge_unoriented <- c()

  for (i in 0:len_p) {
    a <- p_extended[2 * i + 1]
    b <- p_extended[2 * i + 2]

    # Black edge if not consecutive in genome permutation
    if (abs(a - b) != 1) {
      edge_list <- c(edge_list, a, b)
      edge_colors <- c(edge_colors, "black")
      edge_unoriented <- c(edge_unoriented, 1)
    }

    # Gray edge if not adjacent in the list
    i_init <- match(2 * i + 1, p_extended)
    i_end  <- match(2 * i + 2, p_extended)

    if (abs(i_end - i_init) != 1) {
      edge_list <- c(edge_list, 2 * i + 1, 2 * i + 2)
      edge_colors <- c(edge_colors, "gray")
      edge_unoriented <- c(edge_unoriented, (i_init + i_end) %% 2)
    }
  }

  # Build graph and assign attributes
  g <- make_empty_graph(n = len_pe, directed = FALSE)
  g <- add_edges(g, edge_list)
  E(g)$color <- edge_colors
  E(g)$unoriented <- edge_unoriented

  return(g)
}



# Check if two edges interleaves

is_interleaving <- function(p_extended, graph_1, graph_2, edge_1, edge_2) {

  # Get positions of the two ends of each edge in the extended permutation
  e1_nodes <- ends(graph_1, edge_1)
  e2_nodes <- ends(graph_2, edge_2)

  e1_pos <- sort(match(e1_nodes, p_extended))
  e2_pos <- sort(match(e2_nodes, p_extended))

  # Check for interleaving: overlap without containment
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

  # Map components to vertices
  component_vertices <- lapply(comp_ids, function(i) which(comps$membership == i))
  component_subgraphs <- lapply(component_vertices, function(vids) induced_subgraph(g, vids))

  unoriented_flags <- sapply(component_subgraphs, function(cycle) all(E(cycle)$unoriented == 1))
  colors <- ifelse(unoriented_flags, "lightblue", "lightsalmon")
  min_extents <- sapply(component_vertices, function(verts) min(match(verts, p_extended)))
  max_extents <- sapply(component_vertices, function(verts) max(match(verts, p_extended)))

  # Build h with attributes
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
      V(cycle_i)$name <- vertices_i  # Set names for is_interleaving

      gray_i <- E(cycle_i)[E(cycle_i)$color == "gray"]

      for (j in (i + 1):length(h)) {
        vertices_j <- component_vertices[[j]]
        cycle_j <- component_subgraphs[[j]]
        V(cycle_j)$name <- vertices_j  # Same

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

hurdles_count <- function(g, query_sequence_unsig){

  h <- components_graph(g, query_sequence_unsig)
  column_names <- c("membership", "min_extent", "max_extent", "contained",
                    "contains", "hurdle", "superhurdle")

  info <- data.frame(matrix(ncol = length(column_names), nrow = 0))
  colnames(info) <- column_names

  comps <- components(h)

  # Calculate the extent for each unoriented component of h
  line <- 0
  for (i in 1:max(comps$membership)) {
    comp_vertices <- which(comps$membership == i)
    if (all(V(h)$unoriented[comp_vertices])) {
      line <- line + 1
      info[line, "membership"] <- i
      info[line, c("min_extent", "max_extent")] <- c(min(V(h)[comp_vertices]$min_extent),
                                                     max(V(h)[comp_vertices]$max_extent))
      info[line, c("contained", "contains", "hurdle", "superhurdle")] <- c(0, 0, 0, 0)
    }
  }

  if (line == 0) return(info)
  if (line == 1) {
    info[1, "hurdle"] <- 1
    return(info)
  }

  # Check pairwise containment between components
  for (i in 1:(line - 1)) {
    for (j in (i + 1):line) {
      min_i <- info[i, "min_extent"]
      max_i <- info[i, "max_extent"]
      min_j <- info[j, "min_extent"]
      max_j <- info[j, "max_extent"]

      if ((min_i < min_j) & (max_j < max_i)) {
        info[i, "contains"] <- 1
        info[j, "contained"] <- 1
      }
      if ((min_j < min_i) & (max_i < max_j)) {
        info[j, "contains"] <- 1
        info[i, "contained"] <- 1
      }
    }
  }

  # Determine which are hurdles (contain but not contained, or contained but not contain)
  for (i in 1:line) {
    if (info[i, "contains"] == 0) {
      info[i, "hurdle"] <- 1
    }
    if (info[i, "contained"] == 0 && info[i, "contains"] == 1) {
      info[i, "hurdle"] <- 2  # Possible greatest hurdle; analyze later
    }
  }

  if (any(info$hurdle == 2)) {
    possible_greatest <- which(info$hurdle == 2)
    min_hurdles <- min(info[info$hurdle == 1, "min_extent"])
    max_hurdles <- max(info[info$hurdle == 1, "max_extent"])

    for (i in possible_greatest) {
      if (info[i, "min_extent"] > min_hurdles || info[i, "max_extent"] < max_hurdles) {
        info[i, "hurdle"] <- 0
      }
    }

    if (any(info$hurdle == 2)) {
      idx_row <- which(info$hurdle == 2)
      comp_id <- info[idx_row, "membership"]
      comp_vertices <- which(comps$membership == comp_id)
      cycle_indices <- V(h)$cycle_index[comp_vertices]
      info_smallest <- which(info$hurdle == 1)

      # Assume it's a hurdle unless proven otherwise
      info[idx_row, "hurdle"] <- 1
      if (length(info_smallest) == 1) return(info)

      for (cycle in cycle_indices) {
        gray_edges <- E(g)[components(g)$membership == cycle & E(g)$color == "gray"]
        for (j in 1:(length(info_smallest) - 1)) {
          for (k in (j + 1):length(info_smallest)) {
            min_j <- info[info_smallest[j], "min_extent"]
            max_j <- info[info_smallest[j], "max_extent"]
            min_k <- info[info_smallest[k], "min_extent"]
            max_k <- info[info_smallest[k], "max_extent"]

            for (edge in gray_edges) {
              edge_extent <- sort(match(ends(g, edge), query_sequence_unsig))
              condition_1 <- (edge_extent[1] < min_j && edge_extent[2] > max_j &&
                                (edge_extent[1] > max_k || edge_extent[2] < min_k))
              condition_2 <- (edge_extent[1] < min_k && edge_extent[2] > max_k &&
                                (edge_extent[1] > max_j || edge_extent[2] < min_j))

              if (condition_1 || condition_2) {
                info[idx_row, "hurdle"] <- 0
                break
              }
            }
            if (info[idx_row, "hurdle"] == 0) break
          }
          if (info[idx_row, "hurdle"] == 0) break
        }
        if (info[idx_row, "hurdle"] == 0) break
      }
    }
  }

  return(info)
}





##### Identifying Superhurdles #####

superhurdles_count <- function(info, g, query_sequence_unsig){

  # If the number of hurdles is even, it is not necessary to calculate superhurdles because of the definition of fortress
  if (sum(info$hurdle)%%2 == 0) {return(info)}

  # If all unoriented component are hurdles, there is no hurdle protecting an unoriented nonhurdle
  if (sum(info$hurdle) == length(info$hurdle)) {return(info)}

  # For all hurdles: remove one by one and check if some non hurdle turned into a hurdle; if it is the case, then the removed hurdle is a superhurdle
  for (ignored in info$membership[info$hurdle == 1]) {

    # Remove the component from g
    comps <- components(g)
    ignore_vertices <- which(comps$membership == ignored)
    g_reduced <- delete_vertices(g, ignore_vertices)

    # Recompute hurdle info without this component
    info2 <- hurdles_count(g_reduced, p_extended)

    # Compare: did any non-hurdle become a hurdle?
    orig_non_hurdles <- subset(info, hurdle == 0 & membership != ignored)
    new_hurdles <- subset(info2, hurdle == 1)

    if (nrow(new_hurdles) == 0 || nrow(orig_non_hurdles) == 0) next

    overlap <- intersect(orig_non_hurdles$membership, new_hurdles$membership)

    if (length(overlap) > 0) {
      info$superhurdle[info$membership == ignored] <- 1
    }
  }

  return(info)
}





##### Identifying Fortress #####

is_fortress <- function(superhurdles){

  if(sum(superhurdles$hurdles)%%2 == 0){
    return(0)
  }
  if(sum(superhurdles$superhurdle) == sum(superhurdles$hurdle)){
    return(1)
  }
  return(0)

  # A permutation is a fortress if the number of hurdles is odd and all hurdles are superhurdles
  # Return 1 if it is fortress
  # Return 0 otherwise

}





##### Inversion Distance #####

inversionDistance <- function(gb_object){

  p <- permutationVector(gb_object)

  p_extended <- extendedPermutation(p)

  bp_graph <- breakpoint_graph(p_extended)
  hurdles <- hurdles_count(bp_graph, p_extended)
  superhurdles <- superhurdles_count(hurdles, bp_graph, p_extended)

  minimal <- bp_count(p_extended) - cycle_count(bp_graph) + sum(hurdles$hurdle) + is_fortress(superhurdles)

  return(minimal)
}
