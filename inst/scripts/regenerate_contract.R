# Regenerate the bundled Giotto contract from the installed GiottoClass.
# Run from the GiottoConnect package root:
#   Rscript inst/scripts/regenerate_contract.R
pkgload::load_all(".", quiet = TRUE)
ct <- GiottoConnect:::.contract_build()
f <- GiottoConnect:::.contract_write(ct, file.path("inst", "contract"))
message("wrote ", f)
print(structure(ct, class = "giottoContract"))
