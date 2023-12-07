library(SpatialExperiment)
library(STexampleData)
library(scuttle)
library(mgcv)

spe <- Visium_humanDLPFC()
spe <- spe[, colData(spe)$in_tissue == 1]
spe <- logNormCounts(spe, log = FALSE)


# Calcualte library size
spe$total <- counts(spe) |> colSums()

# Smaller set of genes
idx <- which(
  rowData(spe)$gene_name %in%
    c("MOBP", "PCP4", "SNAP25",
      "HBB", "IGKC", "NPY")
)
spe <- spe[idx, ]

coords <- spatialCoords(spe)
colnames(coords) <- c("coord_x", "coord_y")
fit.df <- data.frame(
  coords
)

idx_MOBP <- which(rowData(spe)$gene_name %in% "MOBP")

fit.df$tmp_y <- counts(spe)[idx_MOBP,]

tp_mdl <- gam(
  formula = tmp_y~s(coord_x, coord_y, bs='tp'),
  family = poisson, data = fit.df,
  offset = log(spe$total))

tmp_pred <- predict(tp_mdl, type = "terms")
spe$margin_spatialVar <- tmp_pred[, "s(coord_x,coord_y)"] |> exp()


# 3D Plot -----------------------------------------------------------------
steps <- 50
coord_x <- with(fit.df, seq(min(coord_x), max(coord_x), length = steps) )
coord_y <-  with(fit.df, seq(min(coord_y), max(coord_y), length = steps) )
newdat <- expand.grid(coord_x = coord_x, coord_y = coord_y)
surf_df <- matrix(predict(tp_mdl, newdat), steps, steps)

z <- surf_df
nrz <- nrow(z)
ncz <- ncol(z)

nbcol <- 100

library(viridis)
Lab.palette <- colorRampPalette(viridis(nbcol))
color <- Lab.palette(nbcol)
# Compute the z-value at the facet centres
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]

# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)

persp(
  coord_x, coord_y, surf_df,
  col = color[facetcol], theta = 45,
  border = NA,
  box = FALSE
)





