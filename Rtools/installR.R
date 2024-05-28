library(BSgenome)
library(devtools)


PAPDIR <- "/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/reference_genomes"
#### Installed reference for Pfalc. Dd2
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

#### Installed reference for Pfalc. 3D7
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

#### Reference for Pfalc. 3D7 with supplemental sequences
refDir <- file.path(PAPDIR, "..",
                    "malaria", "cowman_lab", "drug_resistance",
                    "iGP3_Paola2024", "supplementedRef52" )
# Check if the package is already installed
if (!requireNamespace("BSgenome.PfalciparumNF54iGP", quietly = TRUE)) {
  # If not installed, try to install it
  install_local(path = refDir,dependencies = FALSE, INSTALL_opts = '--no-lock')
}
#Check
library("BSgenome.PfalciparumNF54iGP",
        character.only = TRUE
)
print(eval(as.name("BSgenome.PfalciparumNF54iGP")))
