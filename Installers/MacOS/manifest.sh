#!/bin/zsh
# Author: J.Gray
# Email
set -e

source ~/.zshrc

~/.pixi/bin/pixi config set default-channels '["conda-forge", "bioconda"]' --global

~/.pixi/bin/pixi init dashboard

cd dashboard

~/.pixi/bin/pixi add "r-base=4.4.3" "r-seurat" "r-devtools" "r-shiny" "r-tidyverse" "r-shinydashboard" "r-shinyjs" "r-shinyjqui" "r-shinyTree" "r-dt" "r-writexl" "r-data.table" "r-nnls" "r-mcmcpack" "r-e1071" "r-locfdr" "r-quadprog" "r-statmod" "r-corpcor" "r-doparallel" "r-ggally" "r-markdown" "r-xfun" "r-promises"
# ~/.pixi/bin/pixi add "bioconductor-singlecellexperiment" "bioconductor-toast" "bioconductor-genomeinfodb"

echo "✅ R version 4.4.3 installed succesfully!"

# ======= Install libgit2 ======
echo "\tInstalling libgit2"
~/.pixi/bin/pixi add libgit2

# ======= Getting App repository =========
echo "2️⃣\tGetting latest dashboard files from repo"
curl -L "https://github.com/jrdalenberg/UnBlender/archive/refs/heads/remove_unused_packages.zip" -o main.zip

## Extracting files
echo "3️⃣\tExtracting dashboard!"
unzip -X main.zip
rm main.zip

## Removing reminants
mv UnBlender-remove_unused_packages/Unblender-app .
yes | rm -rf UnBlender-remove_unused_packages

# ======= Installing CRAN packages =======
echo "4️⃣\tInstalling required CRAN packages"
# ~/.pixi/bin/pixi run Rscript ../installR_scripts/install_CRAN_packages.R

# ======= Installing BioConductor packages =======
echo "5️⃣\tInstalling required BioConductor packages"
~/.pixi/bin/pixi run Rscript ../installR_scripts/install_BIOC_packages.R

# ======= Installing Github packages =======
echo "5️⃣\tInstalling required Github packages"
~/.pixi/bin/pixi run Rscript ../installR_scripts/install_REPO_packages.R

# ======= Installing runtime documents =====
curl -L "https://ironman.tenwisedev.nl/public/so_downs_5000.rds" -o Unblender-app/app/datafiles/so_downs_5000.rds
curl -L "https://ironman.tenwisedev.nl/public/example_gse76225.csv" -o Unblender-app/app/datafiles/example_gse76225.csv

echo
echo
echo "✅ DASHBOARD INSTALLED SUCCESFULLY!"
