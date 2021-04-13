# newick_parser.R
# Phylogenetic tree reader of newick format file.
# Assumptions: (1) there are no quotation marks (labels),
#              (2) all edges have length
#              (3) file is a single line (no "\n")
# Author: Juan J. Garcia Mesa

suppressMessages(library(stringr))     # string operations

read_newick = function(file) {

  tree <- scan(file, what = "", sep = "\n", quiet = TRUE)

  if (identical(tree, character(0))) {
    warning("Empty character string.")
    return(NULL)
  }

  tree <- gsub("[ \t]", "", tree)   # remove tabs

  tree.length = nchar(tree)
  tree.list = strsplit(tree,"")[[1]]

  ninode = str_count(tree,"\\(")
  nleaf = str_count(tree,",")+1

  edges = NULL    # list of edges: from,to,length
  inode.label = leaf.label = c()
  inode = seq(nleaf+1,nleaf+ninode)

  i = 1
  inode.pos = 0
  current_leaf = 1

  while(i <= tree.length) {

    if(identical(tree.list[i],'(')) { # new inode
      inode.pos = inode.pos + 1
      i = i + 1

    } else if(identical(tree.list[i],')')) { # add inode info and go up

      i = i + 1

      if(identical(tree.list[i],';')) {
        inode.label = c(inode.label,NA)
        break
      } else if(grepl("[A-Za-z]",tree.list[i])) { # read label name if present
        inode.label = c(inode.label, str_extract(substr(tree,i,tree.length),"[A-Za-z0-9_]+"))
        i = i + str_length(inode.label[length(inode.label)]) + 1

      } else {  # insert empty label
        inode.label = c(inode.label, NA)
         i = i + 1
      }

      # read edge length
      len = str_extract(substr(tree,i,tree.length),"[0-9]+(.[0-9]+)?")
      edges = rbind(edges,c(inode[inode.pos-1],inode[inode.pos],as.numeric(len)))

      i = i + nchar(len)
      if(identical(tree.list[i],","))
        i = i + 1

      # update current inode (go up)
      inode = inode[-inode.pos]
      inode.pos = inode.pos - 1

    } else if(identical(tree.list[i],':') | grepl("[A-Za-z]",tree.list[i])) {
      # read label and update current leaf

      if(grepl("[A-Za-z]",tree.list[i])) { # read label name if present
        leaf.label = c(leaf.label, str_extract(substr(tree,i,tree.length),"[A-Za-z0-9_]+"))
        i = i + str_length(leaf.label[length(leaf.label)]) + 1

      } else {  # insert empty label
        leaf.label = c(leaf.label, NA)
      }

      # read edge length
      len = str_extract(substr(tree,i,tree.length),"[0-9]+(.[0-9]+)?")
      edges = rbind(edges,c(inode[inode.pos],current_leaf,as.numeric(len)))

      i = i + nchar(len)
      if(identical(tree.list[i],','))
        i = i + 1

      # update current leaf
      current_leaf = current_leaf + 1

    } else if(identical(tree.list[i],';')) {  # end of tree

      if(i != length(tree.list)) {
        warning("End of tree char ';' is not at end of input string")
        return(NULL)
      }

    } else {
      warning("Character unidentified at position ",i," '", tree.list[i],"'")
      return(NULL)
    }
  }

  # order of alignment of sequences (clades)
  cld = clades(tree,edges,ninode,nleaf)

  # substitute leaf ID with leaf label
  for(i in 1:length(cld)) {
    for(j in 1:length(cld[[i]])) {
      cld[[i]][j] = leaf.label[as.numeric(cld[[i]][j])]
    }
  }

  # return alignment order
  return(cld)
}

# create a list of clades (only containing leafs) that will determine the order
#   in which the sequences are aligned
clades = function(tree, edges, ninode, nleaf) {
  clade.list = list()
  inode = (nleaf+1):(nleaf+ninode)
  for(i in ninode:1) {
    cld = edges[which(edges[,1] == inode[i] & edges[,2] <= nleaf),2]
    if(length(cld) > 0)
      clade.list[[length(clade.list)+1]] = cld
  }
  return(clade.list)
}

if(!interactive()) {
	ARGS = commandArgs(trailingOnly = TRUE)
	read_newick(ARGS[1])
}
