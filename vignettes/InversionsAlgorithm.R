#install.packages("igraph")
library(igraph)

##### Manipulating and creating sequences #####

#transforming a signed permutation into an unsigned permutation to construct the breakpoint graph
sig2unsig <- function(p){
  p_unsig <- sapply(p, function(i) {
    if (i > 0){ c( 2*i,     2*i + 1) }
    else      { c(-2*i + 1, -2*i   ) }
  })
  c(1, p_unsig, 2*length(p) + 2)
}

#perform an inversion
inversion <- function(seq, i, j){
  if  (i%%2 != 0 | j%%2 != 1)               stop("Invalid operation.
                                                      Insert an even starting and an odd ending")
  if (i > length(seq) | j > length(seq))    stop("Out of range index")
  c(seq[1:(i-1)], seq[j:i], seq[-(1:j)])
}
sample_ <- function(x, size) {
  if (length(x) == 1) return(x)
  sample(x, size)
}

#simulate a number of random inversions given by "simulations" in an unsigned sequence s
simulate_inversion <- function(s, simulations){
  n <- length(s)/2 - 2 #length of original sequence
  positions <- matrix(NA, nrow=simulations, ncol=2)

  for (k in 1:simulations) {
    i <- sample(1:n, 1)
    j <- 2*sample_(i:n, 1) + 1
    i <- 2*i
    s <- inversion(s, i, j)  #perform the inversion
    positions[k, ] <- c(i, j)  #store the positions
  }
  list("positions" = positions, "sequence" = s)
}

##### Breakpoint Graph #####

breakpoint_graph <- function(query_sequence_unsig){

  len_qsu <- length(query_sequence_unsig)
  len_qs <- len_qsu/2 -1

  #create an empty graph
  g <- graph.empty(n = len_qsu, directed = FALSE)

  for(i in 0:len_qs) {
    #add black edges between seq[2i+2] and seq[2i+1] if they are not consecutive
    if(abs((query_sequence_unsig[2*i + 2] - query_sequence_unsig[2*i +1])) != 1) {
      g <- add_edges(g, c(query_sequence_unsig[2*i + 2],
                          query_sequence_unsig[2*i +1]),
                     color="black",
                     unoriented=1)}

    #add gray edges between 2i+2 and 2i+1 if they are not adjacent
    i_init <- match((2*i + 1), query_sequence_unsig)
    i_end <- match((2*i + 2), query_sequence_unsig)

    if ((abs(i_end - i_init)) != 1){
      g <- add_edges(g, c(2*i+1,
                          2*i+2),
                     color="gray",
                     unoriented = (i_init+i_end)%%2) #1 if it is unoriented
    }
  }
  return(g)
}


#check if two edges interleaves
is_interleaving <- function(query_sequence_unsig, graph_1, graph_2, edge_1, edge_2) {
  a <- sort(match(c(ends(graph_1, edge_1)[1], ends(graph_1, edge_1)[2]), query_sequence_unsig))
  b <- sort(match(c(ends(graph_2, edge_2)[1], ends(graph_2, edge_2)[2]), query_sequence_unsig))
  return((a[1] < b[1] && b[1] < a[2] && a[2] < b[2]) || (b[1] < a[1] && a[1] < b[2] && b[2] < a[2]))
}

#count the number of breakpoints (including artificial extremities)
bp_count <- function(seq){
  sum(abs(diff(seq)) != 1 )}

#count the number of cycles
cycle_count <- function(g){
  sum(components(g)$csize > 1)
}

##### Components Graph #####

#connected components graph given the breakpoint graph
components_graph <- function(g, query_sequence_unsig){

  cycles <- components(g)

  h <- graph.empty(directed = FALSE)

  #creating nodes of h
  for (i in which(cycles$csize > 1)) {
    #vertices belonging to the current component
    vertices <- which(cycles$membership == i)
    cycle <- induced_subgraph(g, vids = vertices)
    h <- add_vertices(h, 1,
                      cycle_index=i,
                      unoriented=all(E(cycle)$unoriented == 1), #if it is true, the node is unoriented
                      color=ifelse(all(E(cycle)$unoriented == 1), "lightblue", "lightsalmon"),
                      min_extent=min(match(vertices, query_sequence_unsig)),
                      max_extent=max(match(vertices, query_sequence_unsig)))
  }

  #creating edges of h
  if (length(h) > 1){
    for (i in 1:(length(h)-1)) {

      vertices_i <- which(cycles$membership == V(h)$cycle_index[i])
      cycle_i <- induced_subgraph(g, vids = vertices_i)
      V(cycle_i)$name <- vertices_i

      for (j in (i+1):length(h)) {
        vertices_j <- which(cycles$membership == V(h)$cycle_index[j])
        cycle_j <- induced_subgraph(g, vids = vertices_j)
        V(cycle_j)$name <- vertices_j

        #compare all pairs of gray edges of each cycle
        for (k in 1:length(E(cycle_i)[E(cycle_i)$color=="gray"])) {
          for (l in 1:length(E(cycle_j)[E(cycle_j)$color=="gray"])){
            if (is_interleaving(query_sequence_unsig, cycle_i, cycle_j,
                                E(cycle_i)[E(cycle_i)$color=="gray"][k],
                                E(cycle_j)[E(cycle_j)$color=="gray"][l])){
              h <- add_edges(h, c(i, j))
            }
          }
        }
      }
    }#end of creating edges of h
  }
  return(h)
} #end of graph construction

##### Identifying hurdles #####

hurdles_count <- function(g, query_sequence_unsig){

  h <- components_graph(g, query_sequence_unsig)
  column_names <- c("membership", "min_extent", "max_extent", "contained",
                    "contains", "hurdle", "superhurdle")

  info <- data.frame(matrix(ncol = length(column_names), nrow = 0))
  colnames(info) <- column_names

  comps <- components(h)

  #calculate the extent for each unoriented component of h:

  line <- 0

  for (i in 1:max(comps$membership)) {

    comp_vertices <- which(comps$membership == i)

    if (all(V(h)$unoriented[comp_vertices])) {
      info[(line + 1), "membership"] <- i
      line <- line + 1

      u_min <- min(V(h)[comp_vertices]$min_extent)
      u_max <- max(V(h)[comp_vertices]$max_extent)

      info[line, c("min_extent", "max_extent")] <- c(u_min, u_max)
      info[line, c("contained", "contains", "hurdle", "superhurdle")] <- c(0, 0, 0, 0)
    }
  }

  #check number of unoriented components

  if (line == 0) {
    return(info)}

  if (line == 1) {
    info["hurdle"][1] <- 1
    return(info)}

  for (i in 1:(line-1)) {
    for (j in 2:line) {

      min_i <- info[i, "min_extent"]
      max_i <- info[i, "max_extent"]

      min_j <- info[j, "min_extent"]
      max_j <- info[j, "max_extent"]

      if ((min_i < min_j) & (min_j < max_j) & (max_j < max_i)) {
        info[i, "contains"] <- 1
        info[j, "contained"] <- 1
      }

      if ((min_j < min_i) & (min_i < max_i) & (max_i< max_j)) {
        info[j, "contains"] <- 1
        info[i, "contained"] <- 1
      }
    }
  }

  #hurdles never are in "the middle". They need to contain and not being contained,
  # or being contained and not contain

  for (i in 1:line) {
    if ((info[i, "contained"] == 0 || info[i, "contained"] == 1) && info[i, "contains"] == 0){
      info[i, "hurdle"] <- 1
    }
    if (info[i, "contained"] == 0 && info[i, "contains"] == 1){
      info[i, "hurdle"] <- 2 #we need to analyze this case
    }
  }

  if (any(info$hurdle == 2)) {

    #first, check if it covers all other hurdles, otherwise it is not a (greatest) hurdle

    possible_greatest <- which(info$hurdle == 2)

    min_hurdles <-  min(info[info$hurdle == 1, ]["min_extent"])
    max_hurdles <-  max(info[info$hurdle == 1, ]["max_extent"])

    for (i in possible_greatest){

      if (info[i, "min_extent"] > min_hurdles || info[i, "max_extent"] < max_hurdles){
        info[i, "hurdle"] <- 0
      }
    }

    #if there is still someone where hurdle == 2 (it is unique), check the steps below

    if (any(info$hurdle == 2)) {

      #attribute cycle_index from vertices in h allows to access the correspondent cycle in g

      greatest <- which(info$hurdle == 2)
      greatest <- info[greatest, "membership"]

      comp_vertices <- which(comps$membership == greatest)

      idx = V(h)$cycle_index[comp_vertices] #get the memberships of the cycles which compose that component

      info_smallest <- which(info$hurdle == 1)

      #we start assuming it is a hurdle; if we find it dividing two hurdles, we put a zero.

      info[greatest, "hurdle"] <- 1

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

              min_j <- info[j, "min_extent"]  # extent(U')[1]
              max_j <- info[j, "max_extent"]  # extent(U')[2]

              min_k <- info[k, "min_extent"]  # extent(U")[1]
              max_k <- info[k, "max_extent"]  # extent(U")[2]

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
                  info[greatest, "hurdle"] <- 0
                  greatest_check <- TRUE
                }
              }
            }
          }
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

##### Performing the rearrangement by inversions #####


inversions_rearrangement <- function(seq, unsigned=TRUE){

  if(!(unsigned))
  {seq <- sig2unsig(seq)}

  len_qsu <- length(seq)
  len_qs <- len_qsu/2 -1

  bp_graph <- breakpoint_graph(seq)
  hurdles <- hurdles_count(bp_graph, seq)
  superhurdles <- superhurdles_count(hurdles, bp_graph, seq)

  inversions <- matrix(NA, nrow=0, ncol=2)

  loops <- len_qs

  for (loop in 1:loops){

    if (!(all(diff(seq)==1))) {

      for (i in 1:len_qs){
        for (j in i:len_qs){
          k <- 2*i
          l <- 2*j+1

          new_seq <- inversion(seq, k, l)

          if (!(all(diff(new_seq)==1))) {

            new_graph <- breakpoint_graph(new_seq)
            new_hurdles <- hurdles_count(new_graph, new_seq)
            new_superhurdles <- superhurdles_count(new_hurdles, new_graph, new_seq)

            delta = bp_count(new_seq) - cycle_count(new_graph) + sum(new_hurdles$hurdle) + is_fortress(new_superhurdles)                    - (bp_count(seq) - cycle_count(bp_graph) + sum(hurdles$hurdle) + is_fortress(superhurdles))

            if (delta == -1){

              seq <- new_seq
              bp_graph <- new_graph
              hurdles <- new_hurdles
              superhurdles <- new_superhurdles
              inversions <- rbind(inversions, c(k, l))
              #print(c(k, l))
              #print(seq)
            }

          }
          else {
            inversions <- rbind(inversions, c(k, l))
            seq <- new_seq}
        }
      }
    }
  }
  return(inversions)
}
