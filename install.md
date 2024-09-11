# Installation
Here we provide instructions for installation on various systems.

## Ubuntu 14.04 (64 bit), Ubuntu 16.04 (64 bit)
1. Install R
   - add the following line to `/etc/apt/source.list`.
     - `deb https://cran.rstudio.com/bin/linux/ubuntu trusty/`
   ```
   sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
   sudo apt-get update
   sudo apt-get install r-base
   ```
2. Install dependencies for `devtools`
   - `sudo apt-get install libssl-dev`
   - `sudo apt-get install libcurl4-openssl-dev`
3. Install R package dependencies
   ```
   devtools
   RColorBrewer
   gplots
   phytools
   ├──ape
   ├──maps
   ├──Rcpp
   doParallel
   rsvd
   progress
   geiger
   knitr
   RcppArmadillo
   weights
   phangorn
   ```
4. Install impute
   - In R, run the following.
     - `install.packages("BiocManager")`
     - `BiocManager::install("impute")`

5. Install from Github
   ```
   library(devtools)
   install_github("nclark-lab/erc")
   ```

## Win 7 (64 bit), Win 10 (64 bit)
1. Install R
   - `https://cran.r-project.org/bin/windows/base/`
2. Install Rtools
   - `https://cran.r-project.org/bin/windows/Rtools/`
   - add `Rtools\bin` to the `Path`
3. Install R package dependcies
   ```
   devtools
   RColorBrewer
   gplots
   phytools
   ├──ape
   ├──maps
   ├──Rcpp
   doParallel
   rsvd
   progress
   geiger
   knitr
   RcppArmadillo
   weights
   phangorn
   ```

4. Install impute
   - In R, run the following.
     - `install.packages("BiocManager")`
     - `BiocManager::install("impute")`

5. Install from Github
  ```
  library(devtools)
  install_github("nclark-lab/RERconverge")
  ```
