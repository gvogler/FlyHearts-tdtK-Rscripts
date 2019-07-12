# Ubuntu
#
# To install the newest version of R run this in the Linux commandline
# sudo apt install apt-transport-https software-properties-common
# sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
# sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

# sudo apt update
# sudo apt-get install r-base


# sudo apt-get install git

# git clone https://github.com/graukatze/R-scripts-for-tdtK-hearts.git

# Packages necessary to build devtools, ggplot2 etc in R
# sudo apt-get install default-jdk
# sudo apt-get install libxml2-dev
# sudo apt-get install libcurl4-openssl-dev 
# sudo apt-get install libssl-dev 
# sudo apt-get install libcairo2-dev
# sudo apt-get install libpq-dev 
# sudo apt install libgdal-dev
# 
# 
# sudo apt-get install mysql-server postgresql libmariadb-client-lgpl-dev
# 
# # for ggplot2
# sudo apt-get install libudunits2-dev
# sudo apt-get install libfftw3-dev
# sudo apt-get install libtiff5-dev
# sudo apt-get install libgfortran4
# bftools from Bioformats
# sudo cp -rf * /usr/local/bin
# 

# Necessary packages to be installed after R updates

# 
# sudo apt-get install libudunits2-dev
# sudo apt-get install libfftw3-dev
# sudo apt-get install libtiff5-dev
# 
# bftools from Bioformats
# sudo cp -rf * /usr/local/bin
# 

# Necessary packages to be installed after R updates
# Windows 10 specific
# To build packages for R in Windows please install the following


install.packages("rJava", type = "source")
## MacOS: sudo R CMD javareconf

install.packages("devtools", dependencies = TRUE)
install.packages("tidyverse", dependencies = TRUE)
install.packages("xlsx", dependencies = TRUE, type="source", INSTALL_opts=c("--no-multiarch"))
install.packages("openxlsx", dependencies = TRUE)
install.packages("reshape2", dependencies = TRUE)
install.packages("curl", dependencies = TRUE)
install.packages("dplyr", dependencies = TRUE)
install.packages("MASS", dependencies = TRUE)
install.packages("plyr", dependencies = TRUE)
install.packages("quantmod", dependencies = TRUE)
install.packages("zoo", dependencies = TRUE)
install.packages("ggplot2", dependencies = TRUE)
install.packages("TTR", dependencies = TRUE)
install.packages("signal", dependencies = TRUE)
install.packages("pracma", dependencies = TRUE)
install.packages("foreach", dependencies = TRUE)
install.packages("foreach", dependencies = TRUE)
install.packages("broom", dependencies = TRUE)
install.packages("gridExtra", dependencies = TRUE)
install.packages("baseline", dependencies = TRUE)
install.packages("feather", dependencies = TRUE)
install.packages("e1071", dependencies = TRUE)
install.packages('matrixStats', dependencies = TRUE)

# For package fifer:
url <- "https://cran.r-project.org/src/contrib/Archive/fifer/fifer_1.1.tar.gz"
pkgFile <- "fifer_1.1.tar.gz"
download.file(url = url, destfile = pkgFile)

# Install dependencies

install.packages(c('party', 'plotrix', 'randomForestSRC', 'Hmisc', 'fields', 'randomForest', 'xtable'))

install.packages(c('randomForest', 'xtable'))

# Install package
install.packages(pkgs=pkgFile, type="source", repos=NULL)

# Delete package tarball
unlink(pkgFile)


library('devtools', 'rJava')


BiocManager::install("BiocStyle")
BiocManager::install("remotes")
BiocManager::install("EBImage")
BiocManager::install("aoles/RBioFormats", dependencies = TRUE, INSTALL_opts=c("--no-multiarch"))
## INSTALL_opts=c("--no-multiarch") # for Windows 64bit and R-64bit only!

## For older version of Bioconductor
# biocLite("remotes")
# biocLite("EBImage")
# biocLite("aoles/RBioFormats", dependencies = TRUE) 
# biocLite("topGO")
# biocLite("ALL")
# biocLite("hgu95av2.db")



