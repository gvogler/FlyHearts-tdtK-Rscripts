# FlyHearts-tdtK-Rscripts
R-script to analyze fluorescent _Drosophila melanogaster_ hearts
https://zenodo.org/badge/192813672.svg
## Installation

### System Requirements
This script collection was tested on Ubuntu Linux (v18), Windows 10 and macOS 10.14. It requires successful installation of the following software:

 * [R](https://cran.cnr.berkeley.edu/) (3.6.0 or higher; install **only** 64-bit version)
 * [RStudio Desktop Free](https://www.rstudio.com/products/rstudio/download/)
 
 * [Fiji](https://fiji.sc/) (not ImageJ-standalone!)
 * Bioformats command line tools ([_bftools_](https://www.openmicroscopy.org/bio-formats/downloads/))
 
 * [Java JDK or JRE](https://www.oracle.com/technetwork/java/javase/overview/index.html) Version 8, 10 or 11 (**12 will not work**)
 * [Rtools](https://cran.r-project.org/bin/windows/Rtools/) - _Windows only (check 'Add to PATH' **AND** use default installation path, C:\Rtools)_


 **NOTE**:
 
 The executables for **Fiji/ImageJ** and **bftools** _need_ to be in the PATH, i.e. accessible system-wide.
 
 To compile R packages on macOS it might be necessary to install an alternative compiler (gcc7+). A solution is provided e.g. here: [rJava on MacOS Sierra 10.12.15: unsupported option fopenmp (Stackexchange)](https://stackoverflow.com/a/51996290/4154930)
 
