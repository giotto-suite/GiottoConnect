# Converting Seurat, AnnData and SpatialExperiment objects to Giotto

GiottoConnect converts Seurat objects, AnnData (`.h5ad`) files and
Bioconductor SpatialExperiment / SpatialFeatureExperiment objects into
Giotto objects with three functions:

- `seuratToGiotto(sobject)`
- `anndataToGiotto("file.h5ad")`
- `spatialExperimentToGiotto(spe)`

All work the same way. They look at the structure of the input, pick the
bundled *mapping* that fits it (for example “Seurat v5 with an FOV
image”, “scanpy/squidpy AnnData” or “SpatialExperimentIO layout”) and
build the Giotto object from that mapping. Add `source = "<directory>"`
to get a Giotto object backed by a GiottoDisk project instead of an
in-memory one.

``` r

library(GiottoConnect)
```

## Seurat to Giotto

As an example we load the small Xenium lung crop from GiottoData with
Seurat’s own Xenium loader, and run the usual Seurat steps.

``` r

library(Seurat)

xenium_dir <- GiottoData::getSpatialDataset("xenium_mini_lung_ffpe")
s <- LoadXenium(xenium_dir, fov = "fov", segmentations = "cell")
s <- NormalizeData(s) |>
    FindVariableFeatures() |>
    ScaleData() |>
    RunPCA(npcs = 20) |>
    FindNeighbors(dims = 1:20) |>
    RunUMAP(dims = 1:20)
```

Convert it:

``` r

g <- seuratToGiotto(s)
g
#> An object of class giotto 
#> >Active spat_unit:  cell 
#> >Active feat_type:  rna 
#> dimensions    : 377, 1531 (features, cells)
#> [SUBCELLULAR INFO]
#> polygons      : cell 
#> features      : rna NegControlProbe UnassignedCodeword NegControlCodeword 
#> [AGGREGATE INFO]
#> expression -----------------------
#>   [cell][rna] raw normalized
#>   [cell][UnassignedCodeword] raw normalized
#>   [cell][NegControlCodeword] raw normalized
#>   [cell][NegControlProbe] raw normalized
#> spatial locations ----------------
#>   [cell] raw
#> dim reduction --------------------
#>   [cell][rna] pca umap
#> nearest neighbor networks --------
#>   [cell][rna] kNN.pca
#> 
#> 
#> Use objHistory() to see steps and params used
```

The Giotto object holds:

- the counts and normalized values;
- the control assays as their own feature types, named as Giotto’s
  Xenium reader names them (`NegControlProbe`, `NegControlCodeword`,
  `UnassignedCodeword`);
- cell and feature metadata;
- cell centroids and segmentation polygons;
- the transcripts;
- PCA (with loadings) and UMAP;
- the kNN network
  [`FindNeighbors()`](https://satijalab.org/seurat/reference/FindNeighbors.html)
  built (as `kNN.pca`).

``` r

GiottoClass::list_feature_info(g)
#>             feat_info
#>                <char>
#> 1:                rna
#> 2:    NegControlProbe
#> 3: UnassignedCodeword
#> 4: NegControlCodeword
```

### On disk

With `source`, the same conversion writes expression, polygons,
transcripts and networks into a GiottoDisk project:

``` r

project <- file.path(tempdir(), "xenium_from_seurat")
g_disk <- seuratToGiotto(s, source = project)
```

## AnnData to Giotto

`.h5ad` files are read directly from disk (with hdf5r); no python is
needed. GiottoConnect ships a small scanpy-style file to try:

``` r

h5ad <- system.file("extdata", "mini_scanpy.h5ad", package = "GiottoConnect")
g2 <- anndataToGiotto(h5ad)
g2
#> An object of class giotto 
#> >Active spat_unit:  cell 
#> >Active feat_type:  rna 
#> dimensions    : 60, 50 (features, cells)
#> [SUBCELLULAR INFO]
#> [AGGREGATE INFO]
#> expression -----------------------
#>   [cell][rna] raw normalized
#> spatial locations ----------------
#>   [cell] raw
#> spatial networks -----------------
#>   [cell] kNN_network
#> dim reduction --------------------
#>   [cell][rna] pca umap
#> nearest neighbor networks --------
#>   [cell][rna] kNN.pca
#> 
#> 
#> Use objHistory() to see steps and params used
```

From a typical scanpy/squidpy file the conversion takes:

| AnnData | Giotto |
|----|----|
| first integer matrix in `X`, `layers/*` or `raw/X` | raw expression |
| `X`, when it is not integer counts | normalized expression |
| `obs` / `var` (or `raw/var` for `raw/X` counts) | cell / feature metadata |
| `obsm/spatial` (or `X_spatial`, or `obs` x/y columns) | spatial locations |
| `obsm/X_*` | dimension reductions (`X_pca` with its loadings) |
| `obsp/distances` + `uns/neighbors` | kNN network (named from the space it was built in) |
| `obsp/spatial_connectivities` (squidpy) | spatial network |
| `var/feature_types` categories | one feature type per category |

For a large file, convert onto disk; the expression matrices are then
streamed from the file instead of being loaded into memory:

``` r

g2_disk <- anndataToGiotto(h5ad, source = file.path(tempdir(), "from_h5ad"))
```

## SpatialExperiment to Giotto

[`spatialExperimentToGiotto()`](https://giotto-suite.github.io/GiottoConnect/reference/spatialExperimentToGiotto.md)
takes the layouts the Bioconductor readers produce. Here the same Xenium
crop is read with SpatialExperimentIO and normalized with scuttle and
scater:

``` r

library(SpatialExperimentIO)
spe <- readXeniumSXE(xenium_dir)
spe <- spe[, Matrix::colSums(SummarizedExperiment::assay(spe, "counts")) > 0]
spe <- scuttle::logNormCounts(spe) |>
    scater::runPCA(ncomponents = 20)

g3 <- spatialExperimentToGiotto(spe)
g3
#> An object of class giotto 
#> >Active spat_unit:  cell 
#> >Active feat_type:  rna 
#> dimensions    : 377, 1530 (features, cells)
#> [SUBCELLULAR INFO]
#> [AGGREGATE INFO]
#> expression -----------------------
#>   [cell][rna] raw normalized
#>   [cell][NegControlProbe] raw
#>   [cell][UnassignedCodeword] raw
#>   [cell][NegControlCodeword] raw
#> spatial locations ----------------
#>   [cell] raw
#> dim reduction --------------------
#>   [cell][rna] pca
#> 
#> 
#> Use objHistory() to see steps and params used
```

The control probes that SpatialExperimentIO keeps as `altExps` become
their own feature types. PCA loadings and variances (scater’s `rotation`
and `varExplained`) go with the PCA.

The bundled mappings cover:

| Input | Mapping | Converted on top of counts, metadata and reductions |
|----|----|----|
| SpatialExperimentIO ([`readXeniumSXE()`](https://rdrr.io/pkg/SpatialExperimentIO/man/readXeniumSXE.html), [`readCosmxSXE()`](https://rdrr.io/pkg/SpatialExperimentIO/man/readCosmxSXE.html), [`readMerscopeSXE()`](https://rdrr.io/pkg/SpatialExperimentIO/man/readMerscopeSXE.html)) | `spe-basic` | `spatialCoords`; altExps as feature types |
| VisiumIO (`TENxVisium()`, `TENxVisiumHD()`) | `spe-visium` | spot/bin positions converted from pixels to microns |
| SpatialFeatureExperiment (e.g. `readCosMX()`) | `sfe` | cell / nucleus segmentations, transcripts (`txSpots`), spatial graphs; control probes in the main matrix as feature types |

## Giotto to AnnData and back

[`giottoToAnnData()`](https://giotto-suite.github.io/GiottoConnect/reference/giottoToAnnData.md)
writes a spatial unit and feature type of a Giotto object (the active
ones by default) as an `.h5ad` file in scanpy’s layout, without python:
normalized values in `X`, raw counts in `layers/counts`, reductions in
`obsm`, and nearest-neighbour and spatial networks in `obsp`. Networks
keep their Giotto names, so reading the file back gives the same object
parts.

``` r

out <- file.path(tempdir(), "xenium.h5ad")
giottoToAnnData(g, out)

g_back <- anndataToGiotto(out)
GiottoClass::list_nearest_networks(g_back)
#>    spat_unit feat_type nn_type    name
#>       <char>    <char>  <char>  <char>
#> 1:      cell       rna     kNN kNN.pca
```

In scanpy, use a network by its name:
`sc.tl.leiden(adata, neighbors_key = "kNN.pca")`.

## Giotto to SpatialExperiment and back

[`giottoToSpatialExperiment()`](https://giotto-suite.github.io/GiottoConnect/reference/giottoToSpatialExperiment.md)
converts a spatial unit and feature type to a SpatialExperiment: raw
counts in `counts`, normalized values in `logcounts`, reductions under
scater’s names (`PCA`, `UMAP`) and networks as colPairs. When the object
has polygons, transcripts or spatial networks, it builds a
SpatialFeatureExperiment that also holds them (`cellSeg`, `txSpots`,
spatialGraphs); `sfe = "no"` keeps a plain SpatialExperiment.

``` r

spe <- giottoToSpatialExperiment(g)
spe
#> class: SpatialFeatureExperiment 
#> dim: 377 1531 
#> metadata(1): giottoConnect
#> assays(2): counts logcounts
#> rownames(377): ABCC11 ACE2 ... VWA5A VWF
#> rowData names(8): vf_vst_counts_mean vf_vst_counts_variance ...
#>   var.features var.features.rank
#> colnames(1531): afihnpcm-1 afiiicla-1 ... ogadmnan-1 ogaemdkp-1
#> colData names(11): orig.ident nCount_Xenium ... nFeature_ControlProbe
#>   sample_id
#> reducedDimNames(2): PCA UMAP
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : x y
#> imgData names(0):
#> 
#> unit:
#> Geometries:
#> colGeometries: centroids (POINT), cellSeg (MULTIPOLYGON) 
#> rowGeometries: txSpots (MULTIPOINT) 
#> 
#> Graphs:
#> sample01:

g_back <- spatialExperimentToGiotto(spe)
GiottoClass::list_nearest_networks(g_back)
#>    spat_unit feat_type nn_type    name
#>       <char>    <char>  <char>  <char>
#> 1:      cell       rna     kNN kNN.pca
```

## Checking a conversion

Every converted object records what was done.
[`mappingProvenance()`](https://giotto-suite.github.io/GiottoConnect/reference/mappingProvenance.md)
shows the mapping used, what each parameter was bound to, and what could
not be carried over, or was only estimated:

``` r

prov <- mappingProvenance(g)
prov$mapping
#> [1] "seurat-assay5"
prov$losses[, c("part_id", "kind", "field", "note")]
#>    part_id                      kind                           field
#> 1     <NA>                   dropped assays/Xenium/layers/scale.data
#> 2     <NA>                   dropped               graphs/Xenium_snn
#> 3     locs                 estimated                           units
#> 4 polygons                 estimated                        geometry
#> 5   points                 estimated                           units
#> 6       nn recovered_from_parameters                            name
#> 7       nn recovered_from_parameters                        distance
#>                                                                                                            note
#> 1                                                           recomputable from normalized values; no Giotto slot
#> 2                                  only the kNN graph is mapped; sNN can be rebuilt with createNearestNetwork()
#> 3                                                                 FOV coordinates carry no unit; assumed micron
#> 4                                       invalid source rings are repaired (make_valid); their vertices can move
#> 5                                                                 FOV coordinates carry no unit; assumed micron
#> 6 Xenium_nn built by FindNeighbors.Xenium.pca on 'pca'; k = 19 neighbours (k.param 20 includes the cell itself)
#> 7            recomputed in pca dims 1:20 (euclidean); weight = 1 / (1 + distance) replaces the binary adjacency
```

[`readbackCheck()`](https://giotto-suite.github.io/GiottoConnect/reference/readbackCheck.md)
reads the source again and compares it with the Giotto object, part by
part. `MATCH` means identical, `DECLARED` a difference the mapping
announces (for example a unit that had to be assumed). `UNDECLARED`
would mean a real problem:

``` r

rb <- readbackCheck(s, g)
table(rb$status)
#> 
#> DECLARED    MATCH 
#>        1       60
rb[rb$status != "MATCH", ]
#> <giottoReadback> PASS: 0 MATCH, 1 DECLARED, 0 UNDECLARED
#>  status   part_id  aspect   detail                                    
#>  DECLARED polygons geometry 2 of 1531 polygons with a different extent
```

How a mapping is chosen

[`findMapping()`](https://giotto-suite.github.io/GiottoConnect/reference/findMapping.md)
profiles the object and lists every bundled mapping with its status:

``` r

f <- findMapping(profileObject(s), ecosystem = "seurat")
f$candidates
#>                 id version status  coverage specificity
#> 1     seurat-assay   0.1.0   none 0.1111111           1
#> 2    seurat-assay5   0.1.0  match 1.0000000           1
#> 3  seurat-slideseq   0.1.0   none 0.2000000           1
#> 4 seurat-visium-v1   0.1.0   none 0.2500000           1
#> 5    seurat-visium   0.1.0   none 0.1428571           2
#>                                                                     reason
#> 1          param 'assay' could not be bound (Assays value_of: class Assay)
#> 2                                                                  applies
#> 3 param 'image' could not be bound (Images first_of_class: class SlideSeq)
#> 4 param 'image' could not be bound (Images first_of_class: class VisiumV1)
#> 5 param 'image' could not be bound (Images first_of_class: class VisiumV2)
```

Mappings are chosen by structure (an Assay5 with an FOV image, a Visium
image, a scanpy file), not by package version. When a part is stored
somewhere else than usual (for example coordinates in `meta.data`
instead of an FOV), the mapping’s declared fallback is used, and
`mappingProvenance()$alternatives` says so. A mapping can also be passed
explicitly: `seuratToGiotto(s, mapping = "seurat-visium")`.

Coordinates and units

Giotto uses `sdimx = x` and `sdimy = -y` for image-space coordinates (y
pointing down), so converted locations, polygons and transcripts are
flipped in y.

- **Visium:** spot positions are converted from pixels to microns with
  the 100 µm spot pitch (Visium HD: the bin size).
- **Others:** coordinates are taken as stored. Their unit is recorded as
  *estimated*, and the conversion warns when the coordinate range does
  not look like microns, for example scaled to \[0, 1\].

What is not converted yet

- Several samples or FOVs in one object: only the first FOV/image is
  used.
- Images and Visium scale factors.
- Seurat SCT assays.

Data that is present but not converted is reported as a warning when the
conversion runs, and declared losses are listed in
[`mappingProvenance()`](https://giotto-suite.github.io/GiottoConnect/reference/mappingProvenance.md).
