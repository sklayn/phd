barplot(res$eig[,2], names.arg = 1:nrow(res$eig))
drawn <-
c("3", "18", "1", "8", "9", "14", "7", "13", "12", "16", "15"
)
plot.PCA(res, select = drawn, axes = 1:2, choix = 'ind', invisible = 'quali', title = '', cex = cex)
wilks.p <-
structure(9.8743464127193e-08, .Names = "station")
wilks.p
sample = sample(rownames(res$call$X), length(rownames(res$call$X)))
res$call$X = res$call$X[sample,]
res$ind$coord = res$ind$coord[sample[!sample %in% rownames(res$ind.sup$coord)],]
res$ind.sup$coord = res$ind.sup$coord[sample[sample %in% rownames(res$ind.sup$coord)],]
drawn <-
c("3", "18", "1", "8", "9", "14", "7", "13", "12", "16", "15"
)
hab <-
"station"
plotellipses(res, axes = 1:2, invisible = 'quali', select = drawn, keepvar = hab, title = '', cex = cex)
drawn <-
c("silt_clay", "moisture_content", "TOM", "gravel")
plot.PCA(res, select = drawn, axes = 1:2, choix = 'var', title = '', cex = cex)
drawn <-
c("Kraimorie", "Sozopol", "Akin", "Agalina", "Chukalya")
plot.PCA(res, select = drawn, axes = 1:2, choix = 'ind', invisible = c('ind', 'ind.sup'), title = '', cex = cex)
res.hcpc = HCPC(res, nb.clust = -1, graph = FALSE)
drawn <-
c("3", "18", "1", "8", "9", "14", "7", "13", "12", "16", "15"
)
plot.HCPC(res.hcpc, choice = 'map', draw.tree = FALSE, select = drawn, title = '')
dimdesc(res, axes = 1:0)
res.hcpc$desc.var
