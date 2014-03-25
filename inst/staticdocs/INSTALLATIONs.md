## 1. R requirement and installation

R (http://www.r-project.org) is a language and environment for statistical computing and graphics. We assume R (version 3.0.2) has been installed in your local machine. It can be installed following quick instructions below for different platforms (Windows, Mac, and Linux).

* Quick link for `Windows` and `Mac`: [Download R for Windows](http://www.stats.bris.ac.uk/R/bin/windows/base/R-3.0.2-win.exe) and [Download R for Mac](http://www.stats.bris.ac.uk/R/bin/macosx/R-latest.pkg).

* Below are `shell command lines in Terminal` (for both `Linux` and `Mac`):

Assume you have a `ROOT (sudo)` privilege:
    
    sudo su
    wget http://www.stats.bris.ac.uk/R/src/base/R-3/R-3.0.2.tar.gz
    tar xvfz R-3.0.2.tar.gz
    cd R-3.0.2
    ./configure
    make
    make check
    make install
    R # start R

Assume you do not have a ROOT privilege and want R installation under your home directory (below `/home/hfang` should be replaced with yours):

    wget http://www.stats.bris.ac.uk/R/src/base/R-3/R-3.0.2.tar.gz
    tar xvfz R-3.0.2.tar.gz
    cd R-3.0.2
    ./configure --prefix=/home/hfang/R-3.0.2
    make
    make check
    make install
    /home/hfang/R-3.0.2/bin/R # start R

## 2. Install the packages from CRAN, R-Forge and github

Notes: below are `R command lines (NOT shell command lines in Terminal)`.

First, install the dependent packages:

    if(!require(hexbin)) install.packages("hexbin",repos="http://www.stats.bris.ac.uk/R")
    if(!require(ape)) install.packages("ape",repos="http://www.stats.bris.ac.uk/R")
    if(!require(supraHex)) install.packages("supraHex",repos="http://R-Forge.R-project.org", type="source")
    if(!require(igraph)) install.packages("igraph",repos="http://www.stats.bris.ac.uk/R")
    if(!require(devtools)) install.packages("devtools",repos="http://www.stats.bris.ac.uk/R")
    if(!require(arcdiagram)) install_github("arcdiagram",  username="gastonstat")

Then, install the dnet package:

    if(!require(dnet)) install.packages("dnet",repos="http://R-Forge.R-project.org", type="source")