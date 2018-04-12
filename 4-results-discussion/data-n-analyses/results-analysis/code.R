barplot(res$eig[,2], names.arg = 1:nrow(res$eig))
drawn <-
c("3", "18", "1", "13", "14", "15", "2", "16", "8", "7", "9")
plot.PCA(res, select = drawn, axes = 1:2, choix = 'ind', invisible = 'quali', title = '', cex = cex)
wilks.p <-
structure(1.59036097740613e-09, .Names = "station")
wilks.p
sample = sample(rownames(res$call$X), length(rownames(res$call$X)))
res$call$X = res$call$X[sample,]
res$ind$coord = res$ind$coord[sample[!sample %in% rownames(res$ind.sup$coord)],]
res$ind.sup$coord = res$ind.sup$coord[sample[sample %in% rownames(res$ind.sup$coord)],]
drawn <-
c("3", "18", "1", "13", "14", "15", "2", "16", "8", "7", "9")
hab <-
"station"
plotellipses(res, axes = 1:2, invisible = 'quali', select = drawn, keepvar = hab, title = '', cex = cex)
drawn <-
c("sorting", "sand", "mean_grain_size", "gravel", "moisture_content", 
"silt_clay")
plot.PCA(res, select = drawn, axes = 1:2, choix = 'var', title = '', cex = cex)
drawn <-
c("Sozopol", "Kraimorie", "Agalina", "Akin", "Paraskeva")
plot.PCA(res, select = drawn, axes = 1:2, choix = 'ind', invisible = c('ind', 'ind.sup'), title = '', cex = cex)
res.hcpc = HCPC(res, nb.clust = -1, graph = FALSE)
drawn <-
c("3", "18", "1", "13", "14", "15", "2", "16", "8", "7", "9")
plot.HCPC(res.hcpc, choice = 'map', draw.tree = FALSE, select = drawn, title = '')
dimdesc(res, axes = 1:1)
res.hcpc$desc.var
