barplot(res$eig[,2], names.arg = 1:nrow(res$eig))
drawn <-
c("9", "7", "87", "41", "84", "26", "6", "27", "10", "83", "44", 
"86", "72", "18", "40", "8", "16", "19")
plot.PCA(res, select = drawn, axes = 1:2, choix = 'ind', invisible = 'quali', title = '', cex = cex)
wilks.p <-
structure(4.34559841713137e-05, .Names = "station")
wilks.p
sample = sample(rownames(res$call$X), length(rownames(res$call$X)))
res$call$X = res$call$X[sample,]
res$ind$coord = res$ind$coord[sample[!sample %in% rownames(res$ind.sup$coord)],]
res$ind.sup$coord = res$ind.sup$coord[sample[sample %in% rownames(res$ind.sup$coord)],]
drawn <-
c("9", "7", "87", "41", "84", "26", "6", "27", "10", "83", "44", 
"86", "72", "18", "40", "8", "16", "19")
hab <-
"station"
plotellipses(res, axes = 1:2, invisible = 'quali', select = drawn, keepvar = hab, title = '', cex = cex)
drawn <-
c("shoot_count", "bg_biomass_wet", "ag_biomass_wet")
plot.PCA(res, select = drawn, axes = 1:2, choix = 'var', title = '', cex = cex)
drawn <-
c("Gradina", "Otmanli", "Poda")
plot.PCA(res, select = drawn, axes = 1:2, choix = 'ind', invisible = c('ind', 'ind.sup'), title = '', cex = cex)
res.hcpc = HCPC(res, nb.clust = -1, graph = FALSE)
drawn <-
c("9", "7", "87", "41", "84", "26", "6", "27", "10", "83", "44", 
"86", "72", "18", "40", "8", "16", "19")
plot.HCPC(res.hcpc, choice = 'map', draw.tree = FALSE, select = drawn, title = '')
dimdesc(res, axes = 1:1)
res.hcpc$desc.var
