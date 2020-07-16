#!/bin/bash


# Install some common utilities # ? do i need this
sudo apt-get update
sudo apt-get install gcc make zlib1g-dev libbz2-dev liblzma-dev libssl-dev libcurl4-openssl-dev autoconf g++ pkg-config unzip libxml2-dev libncurses5-dev libncursesw5-dev -y


# Install METAL for GWAS meta-analysis
wget http://csg.sph.umich.edu/abecasis/Metal/download/Linux-metal.tar.gz
tar -xvzf Linux-metal.tar.gz



# Install R
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu bionic-cran35/'
sudo apt-get update
sudo apt-get install r-base -y


# Install the required R packages
sudo su - -c "R -e \"install.packages('tidyverse', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('readr', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('dplyr', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('devtools', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('stringr', repos='http://cran.rstudio.com/')\""
sudo su - -c "R -e \"install.packages('forestplot', repos='http://cran.rstudio.com/')\""


sudo apt-get install libgmp-dev

sudo su - -c "R -e \"devtools::install_github('MRCIEU/TwoSampleMR')\""
sudo su - -c "R -e \"install.packages('MendelianRandomization', repos='http://cran.rstudio.com/')\""
