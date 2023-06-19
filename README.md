# mesh-sig-detect
Time to signal analysis for spontaneous report data on pelvic mesh medical device

# Data source

The data is thanks to [curtis-murray](https://github.com/curtis-murray) at his [MedicalDevicesNLP](https://github.com/curtis-murray/MedicalDevicesNLP) repo 


# Repository file structure

* [dat/](https://github.com/tystan/mesh-sig-detect/tree/main/dat) contains the raw and analysis ready data
* [out/](https://github.com/tystan/mesh-sig-detect/tree/main/out) contains processing and analysis output
* [r/](https://github.com/tystan/mesh-sig-detect/tree/main/r) contains `.R` scripts that are prefixed with numbers for rough ordering of steps in processing/analysis/reporting etc
* The `.qmd` files in the top level directory are the [Quarto](https://quarto.org/) (Rmarkdown/Jupiter Notebook-esque documents) that generate the corresponding summary `*.pdf` documents in the top level directory
