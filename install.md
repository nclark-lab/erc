

# Ubuntu 14.04 (64 bit), Ubuntu 16.04 (64 bit)
1. Install R
   - add the following line to `/etc/apt/source.list`.
     - `deb https://cran.rstudio.com/bin/linux/ubuntu trusty/`
   - `sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9`
   ```
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

4. Install from Github
   ```
   library(devtools)
   install_github("nclark-lab/erc")
   ```
