# GiottoConnect

**GiottoConnect** converts objects from other single-cell and spatial ecosystems into Giotto objects.

## Converting Seurat, AnnData and SpatialExperiment objects

```r
library(GiottoConnect)

g <- seuratToGiotto(seurat_object)                 # in memory
g <- anndataToGiotto("sample.h5ad", source = "my_project/")  # GiottoDisk-backed
g <- spatialExperimentToGiotto(spe)                # SpatialExperiment / SFE

giottoToAnnData(g, "sample.h5ad")                   # and back: AnnData
s <- giottoToSeurat(g)                              #           Seurat

mappingProvenance(g)   # mapping used, parameters, losses
```

Each converter profiles the input, picks the bundled mapping that fits its
structure and builds the Giotto object from it. `.h5ad` files are read directly
(hdf5r); no python is needed. See the *Converting Seurat, AnnData and SpatialExperiment objects to Giotto*
article.
