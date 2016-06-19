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
Clone the repository to your local machine and run setup.py installer:

    $: git clone https://github.com/biologyguy/BuddySuite.git
    $: cd BuddySuite
    $: python3 setup.py install
    
This will ensure you have downloaded the necessary python dependencies and adds links to each program into your path.
 For example, the following command should return something similar (the exact path will be system specific):
 
    $: which seqbuddy
    >>> /usr/local/anaconda/bin/seqbuddy

You may choose to also create a set of short-form symbolic links to each program, which can be very convenient if
 using the programs on a regular basis. First check to make sure your short-form command doesn't already exist; the
 following should not return anything:
 
    $: which sb
    >>>

Assuming no conflict, create a new symbolic link (note that the following can be system specific, your installation
 could be in a different directory):

    $: ln -s /usr/local/anaconda/bin/seqbuddy /usr/local/anaconda/bin/sb
    $: which sb
    >>> /usr/local/anaconda/bin/sb

A lighter-weight alternative to creating symbolic links is to create an alias for the command in your shell's rc file (e.g. ~/.bahsrc).
 Use the following example to decide what to add to your rc file:

    # BuddySuite aliases
    alias alb='alignbuddy'
    alias db='databasebuddy'
    alias pb='phylobuddy'
    alias sb='seqbuddy'

Then `which` would return something along the lines of:

    $: which sb
    >>> alias sb='seqbuddy'
    >>>         /usr/local/anaconda/bin/sb

All of the [examples in the wiki](https://github.com/biologyguy/buddysuite/wiki) use the following short forms:

*Tool* | *Short-form*
---------- | -------- 
AlignBuddy | alb
DatabaseBuddy | db
PhyloBuddy | pb
SeqBuddy | sb

## Dependencies
This project has been written in Python3 and is not backwards compatible with Python2. If Python3 is not currently
 installed on your system, I highly recommend using the free [Anaconda manager](http://continuum.io/downloads#py34)
 from Continuum Analytics (if you experience any difficulty, 
 [click here](https://github.com/biologyguy/BuddySuite/wiki/anaconda)). Alternatively, the software can be downloaded 
 directly from the [Python Software Foundation](https://www.python.org/downloads/).

AlignBuddy and PhyloBuddy can be used to launch a number of third party alignment and tree building programs, but
 installation of these optional programs is up to you. For example, if you wish to use PhyloBuddy to build a 
 phylogenetic tree with RAxML, you will first need to get RAxML into your system PATH. 

The SeqBuddy blast, bl2seq, and purge functions require access to the blastp, blastn, and blastdbcmd binaries from the
 [NCBI C++ toolkit](http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/). If not already in your PATH, SeqBuddy.py will
 attempt to download the binaries if any BLAST dependant functions are called.
 
See the [Dependencies](https://github.com/biologyguy/BuddySuite/wiki/Dependencies) page for a full list of all
 third-party packages and software that BuddySuite requires or wraps.
 
## Getting started
Once installed, you can access the modules from the command line using their full names:

    $: seqbuddy -h

Or the short-form shortcuts if you create them:

    $: sb -h

For a detailed breakdown of the tools available within each module, check out the
 [BuddySuite wiki](https://github.com/biologyguy/BuddySuite/wiki).

## Development version installation
The easiest way to get the development environment up and running is to
 [fork](https://help.github.com/articles/fork-a-repo/) the repository and then clone from your github account:

    $: git clone https://github.com/<YOUR USER ID>/BuddySuite.git

Then move into the repo, switch to the 'develop' branch, and run setup.py (to get dependencies):
    
    $: cd BuddySuite
    $: git checkout develop
    $: python3 setup.py install

Now set the symbolic links to the repo version of each tool, e.g.:

    $: ln -s /<your>/<path>/<to>/BuddySuite/buddysuite/SeqBuddy.py /usr/local/bin/sb
    $: which sb
    >>> /usr/local/bin/sb

All of the individual Buddy toolkits are located in the 'buddysuite' directory. The 
 ['develop' branch](https://github.com/biologyguy/BuddySuite/tree/develop) is where all new features are created
 and tested, so things may be less stable here; it's usually pretty solid though. If you're interested in contributing
 to the project, please ensure you are working from this branch.

See the [developer page](https://github.com/biologyguy/BuddySuite/wiki/Developers) for further information on
 development version dependencies and how to contribute to the project.

## Unit tests
We are striving for high unit test coverage with py.test. There are two ways to run the unit tests, each of which
 should be executed before making a pull request. The first method is faster and will be used more frequently:

    $: cd BuddySuite/buddysuite
    $: bash run_tests.hs

The second method should be run just before submitting a pull request, and uses
 [Docker](https://docs.docker.com/engine/installation/) to build a clean Linux environment with all dependencies. 
 The tests are then run in a container that mimics the environment used by Travis-CI to monitor the state of the BuddySuite
 repository:
    
    $: cd BuddySuite
    $: docker build -t docker-build:latest docker-build
    $: docker run -v $PWD:/home/docker/BuddySuite docker-build:latest

## Citation
There is a very short application note on bioRxiv that can be cited if you use BuddySuite in your work.

[DOI: 10.1101/040675](http://dx.doi.org/10.1101/040675)


## Contact
Any comments you may have would be really appreciated. Please feel free to add issues in the GitHub issue tracker or
 contact me directly at [steve.bond@nih.gov](mailto:steve.bond@nih.gov)
