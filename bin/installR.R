library(BSgenome)
library(devtools)


PAPDIR <- "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/reference_genomes"

refDir <- file.path( 
    PAPDIR, "plasmodium", "PlasmoDB-57_PfalciparumDd2","BSgenome.PfalciparumDd2.PlasmoDB.57")

# Check if the package is already installed
if (!requireNamespace("BSgenome.PfalciparumDd2.PlasmoDB.57", quietly = TRUE)) {
  # If not installed, try to install it
  install_local(path = refDir,dependencies = FALSE, INSTALL_opts = '--no-lock')
}
#Check
library("BSgenome.PfalciparumDd2.PlasmoDB.57",
  character.only = TRUE
)
print(eval(as.name("BSgenome.PfalciparumDd2.PlasmoDB.57")))


refDir <- file.path(PAPDIR,
                    "plasmodium","PlasmoDB-52_Pfalciparum3D7","BSgenome.Pfalciparum3D7.PlasmoDB.52")
# Check if the package is already installed
if (!requireNamespace("BSgenome.Pfalciparum3D7.PlasmoDB.52", quietly = TRUE)) {
  # If not installed, try to install it
  install_local(path = refDir,dependencies = FALSE, INSTALL_opts = '--no-lock')
}
#Check
library("BSgenome.Pfalciparum3D7.PlasmoDB.52",
  character.only = TRUE
)
print(eval(as.name("BSgenome.Pfalciparum3D7.PlasmoDB.52")))


refDir <- file.path(PAPDIR,
                    "plasmodium","vivaxP01","BSgenome.PvivaxP01.PlasmoDB.52")
# Check if the package is already installed
if (!requireNamespace("BSgenome.PvivaxP01.PlasmoDB.52", quietly = TRUE)) {
  # If not installed, try to install it
  install_local(path = refDir,dependencies = FALSE, INSTALL_opts = '--no-lock')
}
#Check
library("BSgenome.PvivaxP01.PlasmoDB.52",
  character.only = TRUE
)
print(eval(as.name("BSgenome.PvivaxP01.PlasmoDB.52")))