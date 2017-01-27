## 1. R requirement

R (http://www.r-project.org) is a language and environment for statistical computing and graphics. We assume R (version 3.1.0 or higher) has been installed in your local machine. The latest version can be installed following instructions below for different platforms (Windows, Mac, and Linux).

* Quick link for `Windows`: [Download R for Windows](http://cran.r-project.org/bin/windows/base).
* Quick link for `Mac`: [Download R for Mac OS X 10.6 (Snow Leopard or higher)](http://cran.r-project.org/bin/macosx).

* Below are `shell command lines in Terminal` (for `Linux`):

Assume you have a `ROOT (sudo)` privilege:
    
    sudo su
    # here enter your password
    wget http://www.stats.bris.ac.uk/R/src/base/R-3/R-3.2.4.tar.gz
    tar xvfz R-3.2.4.tar.gz
    cd R-3.2.4
    ./configure
    make
    make check
    make install
    R # start R

Assume you do not have a ROOT privilege and want R installation under your home directory ($HOME):

    wget http://www.stats.bris.ac.uk/R/src/base/R-3/R-3.2.4.tar.gz
    tar xvfz R-3.2.4.tar.gz
    cd R-3.2.4
    ./configure --prefix=$HOME/R-3.2.4
    make
    make check
    make install
    $HOME/R-3.2.4/bin/R # start R

## 2. Installation of the package

Notes: below are `R command lines (NOT shell command lines in Terminal)`.

First, install dependant/imported/suggested packages:

    source("http://bioconductor.org/biocLite.R")
    biocLite(c("hexbin","ape","supraHex","graph","Rgraphviz","igraph","Biobase","limma","survival","foreach","doMC","devtools"))

Second, install the package `dnet` under [stable release version hosted in CRAN](http://cran.r-project.org/package=dnet):

    install.packages("dnet",repos="http://cran.r-project.org",type="source")

Third (`highly recommended`), update the package `dnet` from [latest development version hosted in GitHub](https://github.com/hfang-bristol/dnet):

    library(devtools)
    if("dnet" %in% rownames(installed.packages())) remove.packages("dnet")
    install_github(c("hfang-bristol/dnet"))

## 3. Workflow of the package

It provides a brief overview of how the package operates and what you expect to get from it.

<a href="javascript:newWin('dnet_workflow.png', 'dnet_workflow.png', '1200', '1200')" title="Click to enlarge"><img style="max-width:95%;border:1px solid #0000FF;box-shadow:5px 5px 2px #C0C0FF;" src='dnet_workflow.png', width="800px" /></a>
