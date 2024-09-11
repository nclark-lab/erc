# Installation
Here we provide instructions for installation on various systems.

## Ubuntu 14.04 (64 bit), Ubuntu 16.04 (64 bit)
1. Install R
   - add the following line to `/etc/apt/sources.list`.
     - `deb https://cran.rstudio.com/bin/linux/ubuntu trusty/`
   - Run the following code in the terminal:
   ```
   sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
   sudo apt-get update
   sudo apt-get install r-base
   ```
2. Install dependcies for `devtools`
   - `sudo apt-get install libssl-dev`
   - `sudo apt-get install libcurl4-openssl-dev`
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

5. Run this command in the terminal to copy the Github repository to your computer.
   ```
   git clone https://github.com/nclark-lab/erc
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
  install_github("nclark-lab/erc")
  ```

## Mac
For compatibility reasons the Mac installation involves some extra dependencies, so be careful to include them all. Also ensure your R is up to date (>4.4.0).

1. Make sure your R version is at least the version under which the binary was compiled. If not, install the latest version of R
   - `https://cran.r-project.org/bin/macosx/`

2. Install the dependent packages as follows:
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
   These packages are unnique to mac:
   ```
   gfortran
   xcodecli
   ```
   - If you have a problem with these, try visiting this site for help:
   - `https://github.com/coatless-mac/macrtools`

3. In the console, navigate to your target folder, and run the following to copy the erc repository to your computer:
   - `git clone https://github.com/nclark-lab/erc`
