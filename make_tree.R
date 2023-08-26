library(ape)
args <- commandArgs(trailingOnly = TRUE)

te <- ape::read.tree(sprintf("%s.treefile", args[1]))

te <- root(te, outgroup = "Lunularia_cruciata", resolve.root = 1)

v = rep("Black", length(te$tip.label))
v[which(te$tip.label == args[1])] = "Red"

pdf(file=sprintf("%s_tree.pdf", args[1]), width = 40, height = 100)
ape::plot.phylo(te, cex = 1.0, tip.color = v)

dev.off()

write.tree(te, sprintf("%s_tree.rooted", args[1]))
