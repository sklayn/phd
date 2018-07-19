barplot(res$eig[,2], names.arg = 1:nrow(res$eig))
drawn <-
c("2", "26", "7", "1", "3", "10", "5", "18", "14", "21", "19", 
"11")
plot.PCA(res, select = drawn, axes = 1:2, choix = 'ind', invisible = 'quali', title = '', cex = cex)
wilks.p <-
structure(8.22559958029619e-06, .Names = "station")
wilks.p
sample = sample(rownames(res$call$X), length(rownames(res$call$X)))
res$call$X = res$call$X[sample,]
res$ind$coord = res$ind$coord[sample[!sample %in% rownames(res$ind.sup$coord)],]
res$ind.sup$coord = res$ind.sup$coord[sample[sample %in% rownames(res$ind.sup$coord)],]
drawn <-
c("2", "26", "7", "1", "3", "10", "5", "18", "14", "21", "19", 
"11")
hab <-
"station"
plotellipses(res, axes = 1:2, invisible = 'quali', select = drawn, keepvar = hab, title = '', cex = cex)
drawn <-
c("Ntotal", "LUSI", "chl_a")
plot.PCA(res, select = drawn, axes = 1:2, choix = 'var', title = '', cex = cex)
drawn <-
c("Poda", "Gradina", "Otmanli")
plot.PCA(res, select = drawn, axes = 1:2, choix = 'ind', invisible = c('ind', 'ind.sup'), title = '', cex = cex)
drawn <-
c("25", "18", "3", "26", "4", "6", "14", "11", "9", "28", "2", 
"16", "24", "8", "12", "19", "10")
plot.PCA(res, select = drawn, axes = 3:4, choix = 'ind', invisible = 'quali', title = '', cex = cex)
wilks.p <-
structure(0.0688114100435798, .Names = "station")
wilks.p
sample = sample(rownames(res$call$X), length(rownames(res$call$X)))
res$call$X = res$call$X[sample,]
res$ind$coord = res$ind$coord[sample[!sample %in% rownames(res$ind.sup$coord)],]
res$ind.sup$coord = res$ind.sup$coord[sample[sample %in% rownames(res$ind.sup$coord)],]
drawn <-
c("25", "18", "3", "26", "4", "6", "14", "11", "9", "28", "2", 
"16", "24", "8", "12", "19", "10")
hab <-
"station"
plotellipses(res, axes = 3:4, invisible = 'quali', select = drawn, keepvar = hab, title = '', cex = cex)
drawn <-
c("NH4", "chl_a", "secchi")
plot.PCA(res, select = drawn, axes = 3:4, choix = 'var', title = '', cex = cex)
drawn <-
c("Ropotamo", "Otmanli")
plot.PCA(res, select = drawn, axes = 3:4, choix = 'ind', invisible = c('ind', 'ind.sup'), title = '', cex = cex)
res.hcpc = HCPC(res, nb.clust = -1, graph = FALSE)
drawn <-
c("2", "26", "7", "1", "3", "10", "5", "18", "14", "21", "19", 
"11")
plot.HCPC(res.hcpc, choice = 'map', draw.tree = FALSE, select = drawn, title = '')
dimdesc(res, axes = 1:3)
res.hcpc$desc.var
