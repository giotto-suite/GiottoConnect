# Follow-ups in other packages (to assess in a separate session)

These were found while benchmarking the GiottoConnect mapping prototype. None is
confirmed as a bug yet. Each entry gives a reproduction, what is established, and the
question that decides between "bug" and "intended behaviour".

## Environment for every reproduction

- **Scratch library:** `R_LIBS=~/Documents/Software/rlibs/suite_dev`, containing:
  - GiottoClass 0.7.31 (`gsource` 6a86287)
  - Giotto 4.3.6 (`gsource` 611c419)
  - GiottoDisk 0.0.0.8 (`dev` f0f431f)
  - GiottoData 0.3.8 (user library)
- **Options:** `options(giotto.has_conda = FALSE, arrow.unsafe_metadata = TRUE)`.
- **Crops:** the GiottoData crops, cached at `~/Library/Caches/org.R-project.R/R/GiottoData/`
  (`getSpatialDataset()`).

---

## 1. GiottoDisk `importXeniumDisk()` drops all-zero genes

**Observed.** The crop `xenium_mini_lung_ffpe` (1,531 cells, 377 gene targets) gives:
- `importXeniumDisk()`: 376 rna features;
- `Seurat::LoadXenium()` on the same folder: 377.

The missing gene is `TCL1A`, which has a total count of 0 in the crop. On the shared genes
the raw counts are identical.

**Reproduce**
```r
library(GiottoDisk)
C <- file.path(GiottoData::giottoDataCache(), "xenium_mini_lung_ffpe")
g <- importXeniumDisk(C, backend = tempfile())$create_gobject()
f_disk <- GiottoClass::featIDs(g, feat_type = "rna")
f_seu  <- rownames(Seurat::LoadXenium(C, fov = "fov")[["Xenium"]])
setdiff(f_seu, f_disk)                     # "TCL1A"
m <- Seurat::LayerData(Seurat::LoadXenium(C, fov = "fov"), layer = "counts")
sum(m["TCL1A", ])                          # 0
```

**Likely cause to check.** The `parquetExprStore` feature ids may be derived from the
stored triplets, so a feature with no non-zero entries leaves no row. The alternative is
that they are taken from `features.tsv` / the h5 `features` table.

**Decide**
- Should a disk-backed expression object keep the full feature panel (as the h5 matrix
  and Seurat do), or only observed features?
- If the panel should be kept: check whether the `exprInput` iterator or the store
  summary drops ids, and whether other importers behave the same (CosMx, Visium HD,
  Atera, Stereo-seq).
- Does any downstream step (feature metadata joins, `featIDs` vs `fDataDT` agreement)
  break when they differ?

---

## 2. GiottoDisk `importVisiumHDDisk()` errors on the crop root

**Observed.** On the crop root the importer warns and then fails:

```r
GiottoDisk::importVisiumHDDisk(file.path(C, "visiumhd_mini_crc"), backend = tempfile(),
    bin = 16, micron = TRUE)
# warnings: "008um_outputs.tar not detected", "square_008um$ not detected", ...
# Error in dir.exists(p) : invalid filename argument
#   (.visiumhd_detect_paths -> .visiumhd_is_segmentation_output -> dir.exists)
```

The same call on `file.path(C, "visiumhd_mini_crc", "binned_outputs")` works. The crop's
own test script (`spatial-datasets/crops/visiumhd_mini_crc/test_visiumhd_mini_crc.R`)
passes `binned_outputs/`.

Notes:
- `bin = 16` was given, but the path detection still looks for the 8 µm outputs.
- The final error comes from `dir.exists()` receiving a non-character (probably `NULL`
  or `NA`) path.

**Decide**
- Is the crop root (the Space Ranger `outs/` equivalent, containing `binned_outputs/`
  and `segmented_outputs/`) a supported input? If yes, path detection should find
  `binned_outputs/square_016um/` for `bin = 16`.
- Independently, `.visiumhd_is_segmentation_output()` should handle a missing path
  without an internal `dir.exists` error, and say which path it expected.

---

## 3. GiottoData: CosMx crop size vs its description

**Observed.** The GiottoData description of `cosmx_mini_colon_wtx` says 9 FOVs and
16,766 cells × 20,378 features. The cached files contain:
- `S0_metadata_file.csv.gz`: 1,752 cells, all in FOV 71;
- `S0_fov_positions_file.csv.gz`: 94 bytes, one row (FOV 71).

`Seurat::LoadNanostring()` therefore loads 1,752 cells.

**Reproduce**
```r
C <- file.path(GiottoData::giottoDataCache(), "cosmx_mini_colon_wtx")
md <- data.table::fread(file.path(C, "S0_metadata_file.csv.gz"), select = c("fov", "cell_ID"))
nrow(md); unique(md$fov)                   # 1752, 71
data.table::fread(file.path(C, "S0_fov_positions_file.csv.gz"))   # 1 row
```

**Decide**
- Which is out of date: the description (in `datasets.json`, the crop README or NEWS), or
  the zip uploaded under the release tag?
- Compare the cached zip's sha256 with `datasets.json`, and with the crop manifest's
  recorded window/selection.
- Note: the crop README and NEWS also said the crop has no transcript file, but the cache
  has a 33 MB `S0_tx_file.csv.gz`.
- If the zip is right, fix the description; if it is a stale upload, re-run the crop script
  (the GiottoData#75 crop code is referenced but not in the repo).

---

## 4. GiottoConnect `giottoToSpatialExperiment()` calls a removed GiottoClass function

**Observed.** It fails with
`could not find function "get_spatial_locations"`. GiottoClass removed the internal
snake_case getters. Before the Seurat tests were ported this was hidden, because the
interop test file failed earlier on `giottoToSeuratV5`.

**Decide.** Port it to `getSpatialLocations()`, or move SpatialExperiment onto the
mapping engine (the next adapter after AnnData).
