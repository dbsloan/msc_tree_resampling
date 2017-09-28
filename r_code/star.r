args <- commandArgs(TRUE)

if (length(args) != 3) stop ("Required command line arugments are: input_tree_file outgroup_name output_name")

library(phybase, warn.conflicts=FALSE, quietly=TRUE)

genetrees<-read.tree.string(file=args[1],format="phylip")

#read.tree.string only extract taxa names from the first tree. Add any additional names from other trees 
for (input_tree in genetrees$tree){
	input_tree_nodes = read.tree.nodes(input_tree)
	for (input_name in input_tree_nodes$names){
		if(!input_name %in% genetrees$names){
			genetrees$names = append(genetrees$names, input_name)
		}
	}
}

#define a species matrix. This assumes that there is one terminal for each species. STAR can accomodate mutliple terminals per species, but the following code would have to be modified to specify that in the matrix.
numtax=length(genetrees$names)
species.structure<-matrix(0,numtax,numtax)
diag(species.structure)<-1

#generate species tree with STAR based on all gene trees. Note this requires specification of a single terminal as an outgroup for rooting the resulting trees. This is specificied by the user on the command line 
full_sptree = star.sptree(genetrees$tree, speciesname=genetrees$names, taxaname=genetrees$names, species.structure=species.structure, outgroup=args[2], method="nj")

#print to output file
sink(args[3], append=TRUE)
cat (c(full_sptree, "\n"), sep="")
