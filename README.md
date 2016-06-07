[![Build Status](https://travis-ci.org/biologyguy/BuddySuite.svg?branch=develop)](https://travis-ci.org/biologyguy/BuddySuite)
<p align="center"><a href="https://github.com/biologyguy/BuddySuite/wiki">
<img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/workshop/images/BuddySuite-logo.gif" /></a></p>
<p align="center">
<a href="https://github.com/biologyguy/BuddySuite/wiki/SeqBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/workshop/images/SeqBuddy-logo.gif" width=20%/></a>
<a href="https://github.com/biologyguy/BuddySuite/wiki/AlignBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/workshop/images/AlignBuddy-logo.gif" width=25%/></a>
<a href="https://github.com/biologyguy/BuddySuite/wiki/DatabaseBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/workshop/images/DBBuddy-logo.gif" width=20%/></a>
<a href="https://github.com/biologyguy/BuddySuite/wiki/PhyloBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/workshop/images/PhyloBuddy-logo.gif" width=25%/></a>
</p>
<p align="center">Do fun stuff with biological data files. Seriously, biological data is fun stuff :)</p>
___
## Description
The BuddySuite modules are designed to be 'one-stop-shop' command line tools for common biological data file
 manipulations.

[SeqBuddy](https://github.com/biologyguy/BuddySuite/wiki/SeqBuddy) is the most mature BuddySuite tool, although
 [AlignBuddy](https://github.com/biologyguy/BuddySuite/wiki/AlignBuddy) and
 [PhyloBuddy](https://github.com/biologyguy/BuddySuite/wiki/PhyloBuddy) are also functional with a more limited number
 of commands. [DatabaseBuddy](https://github.com/biologyguy/BuddySuite/wiki/DatabaseBuddy) is a very different project,
 existing mostly as a 'live shell' for downloading sequences from GenBank, ENSEMBL, and UniProt.

Being pure Python, the BuddySuite should be cross platform. Development and testing have been done on Linux
 and Mac OS X, however, so it is unclear how well the suite will work within Windows.

## Installation 
Installation should be extremely easy on Mac and Linux using the graphical installer (Windows users must install the
 development version, [see below](https://github.com/biologyguy/BuddySuite#development-version-installation)).

[Click here](https://raw.github.com/biologyguy/BuddySuite/master/BuddySuite_installer.py) to download the graphical
 installer and run it from the command line
    
    $: cd /path/to/download/folder
    $: chmod +x BuddySuite_installer.py
    $: ./BuddySuite_installer.py

By default, the installer will create short-form symbolic links for the main tools in your system $PATH ('sb' for
 SeqBuddy, 'alb' for AlignBuddy, 'pb' for PhyloBuddy, and 'db' for DatabaseBuddy), so they can be accessed quickly
 ([examples in the wiki](https://github.com/biologyguy/buddysuite/wiki) use these short forms). The full names of each
 tool will also be added to $PATH. If working outside the context of a graphical OS (on a cluster, for example), the
 installer will run in command-line mode (also accessible with the -cmd flag on graphical systems, if you prefer that).

The BuddySuite installer will only bundle stable release versions of the BuddySuite. If bugs are found they will be
 hot-fixed, but the *expected* behavior will not be changed once the release is finalized. Likewise, new features added
 to the development versions will not become available in the installer until the next release. Versions of each tool or
 the installer can be displayed using the -v flag.

## Dependencies
This project has been written in Python3 and is not backwards compatible with Python2. If Python3 is not currently
 installed on your system, I highly recommend using the free [Anaconda manager](http://continuum.io/downloads#py34)
 from Continuum Analytics (if you experience any difficulty, 
 [click here](https://github.com/biologyguy/BuddySuite/wiki/anaconda)). Alternatively, the software can be downloaded 
 directly from the [Python Software Foundation](https://www.python.org/downloads/).

AlignBuddy and PhyloBuddy can be used to launch a number of third party alignment and tree building programs, but
 installation of these optional programs is up to you. For example, if you wish to use PhyloBuddy to build a 
 phylogenetic tree with RAxML, you will first need to get RAxML into your system PATH. 

All other dependencies come prepackaged with the installer, so you only need to worry about the following if you
 are using the unstable workshop version of BuddySuite.

The SeqBuddy blast, bl2seq, and purge functions require access to the blastp, blastn, and blastdbcmd binaries from the
 [NCBI C++ toolkit](http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/). If not already in your PATH, SeqBuddy.py will
 attempt to download the binaries if any BLAST dependant functions are called. [BioPython](http://biopython.org/) is
 used heavily by the entire suite; any version earlier than 16.6 will cause unit tests to fail. PhyloBuddy requires 
 [DendroPy](https://pythonhosted.org/DendroPy/) and version 3.0 (beta) of the
 [ETE toolkit](http://etetoolkit.org/download/).
 
## Getting started
Once installed, you can access the modules from the command line using their full names:

    $: SeqBuddy -h

Or the shortcuts created by the installer:

    $: sb -h

For a detailed breakdown of the tools available within each module, check out the
 [BuddySuite wiki](https://github.com/biologyguy/BuddySuite/wiki).

## Development version installation
The easiest way to get the development version up and running is to
 [clone/fork](https://help.github.com/articles/fork-a-repo/) the repository.

    $: git clone https://github.com/biologyguy/BuddySuite.git

Then move into the repo and switch to the 'development' branch:
    
    $: cd BuddySuite
    $: git checkout develop

All of the individual Buddy toolkits are located in the 'buddysuite' directory. The 
 ['develop' branch](https://github.com/biologyguy/BuddySuite/tree/develop) is where all new features are created
 and tested, so things may be less stable here; it's usually pretty solid though. If you're interested in contributing
 to the project, please ensure you are working from this branch.

To run the tests, use docker to get all of the dependencies installed:

    $: docker build -t docker-build:latest docker-build

and then run the tests with

    $: docker run -v $PWD:/home/docker/BuddySuite docker-build:latest

See the [developer page](https://github.com/biologyguy/BuddySuite/wiki/Developers) for further information on
 development version dependencies and how to contribute to the project.


## Citation
There is a very short application note on bioRxiv that can be cited if you use BuddySuite in your work.

[DOI: 10.1101/040675](http://dx.doi.org/10.1101/040675)


## Contact
Any comments you may have would be really appreciated. Please feel free to add issues in the GitHub issue tracker or
 contact me directly at [steve.bond@nih.gov](mailto:steve.bond@nih.gov)
