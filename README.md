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
 
If you're new to the command line, or simply want to get a better feel for how BuddySuite works, check out the [Beginners Guide](https://github.com/biologyguy/BuddySuite/wiki/Beginners-Guide).

## Installation 
Clone the repository to your local machine and run setup.py installer (or
 [fork](https://help.github.com/articles/fork-a-repo/) it if planning to develop):

    $: git clone https://github.com/biologyguy/BuddySuite.git
    $: cd BuddySuite
    $: python3 setup.py install clean
    
The project should now be in your Python PATH along with sym-links to the executables. For example:
 
    $: which seqbuddy
    >>> /usr/local/anaconda/bin/seqbuddy

And to get started, simply use the 'help' flag

    $: seqbuddy -h

For a detailed breakdown of the tools available within each module, check out the
 [BuddySuite wiki](https://github.com/biologyguy/BuddySuite/wiki).
 
### Create shortcuts (OPTIONAL)
We also recommend creating a set of short-form symbolic links or aliases to each program, which are very convenient if
 using the programs on a regular basis. Furthermore, all of the 
 [examples in the wiki](https://github.com/biologyguy/buddysuite/wiki) use the following short forms:

*Tool* | *Short-form*
---------- | -------- 
AlignBuddy | alb
DatabaseBuddy | db
PhyloBuddy | pb
SeqBuddy | sb
 
First check to make sure the short-form command doesn't already exist; for example, the
 following should not return anything:
 
    $: which sb
    >>>

If there is a conflict, choose something else that seems reasonable (e.g., 'sbd' or 'sqb'). Now either create aliases
 or symbolic links:

**Aliases**

Copy the following into either the `.bashrc`, `.profile`, or `.bash_profile` (system dependent) file in your home directory:
  
    alias sb="seqbuddy"
    alias alb="alignbuddy"
    alias pb="phylobuddy"
    alias db="databasebuddy"

Then restart your terminal

**Symbolic Links**

For each program, replace the necessary path components in the command below and run from the shell

    $: ln -s /<your>/<path>/<to>/BuddySuite/buddysuite/SeqBuddy.py /usr/local/bin/sb

You should now be able to see the short form commands in your PATH, e.g.,:

    $: which sb
    >>> /usr/local/bin/sb

## Dependencies
This project has been written in Python3 and is not backwards compatible with Python2. If Python3 is not currently
 installed on your system, we highly recommend using the free [Anaconda manager](http://continuum.io/downloads#py34)
 from Continuum Analytics (if you experience any difficulty, 
 [click here](https://github.com/biologyguy/BuddySuite/wiki/anaconda)). Alternatively, the software can be downloaded 
 directly from the [Python Software Foundation](https://www.python.org/downloads/).

AlignBuddy and PhyloBuddy can be used to launch a number of third party alignment and tree building programs, but
 installation of these optional programs is up to you. For example, if you wish to use PhyloBuddy to build a 
 phylogenetic tree with RAxML, you will first need to get RAxML into your system PATH. 

The SeqBuddy blast, bl2seq, and purge functions require access to the blastp, blastn, and blastdbcmd binaries from the
 [NCBI C++ toolkit](http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/).
 
See the [Dependencies](https://github.com/biologyguy/BuddySuite/wiki/Dependencies) page for a full list of all
 third-party packages and software that BuddySuite requires or wraps.
 
## Developers
All of the individual Buddy toolkits are located in the 'buddysuite' directory and the 
 ['develop' branch](https://github.com/biologyguy/BuddySuite/tree/develop) is where all new features have been
 implemented. If you're interested in contributing, please refer to the
 [developer page](https://github.com/biologyguy/BuddySuite/wiki/Developers) for further information on dependencies
 and instructions.

## Unit tests
We are striving for high unit test coverage with py.test. There are two ways to run the unit tests, each of which
 should be executed before making a pull request. The first method is faster and will be used more frequently:

    $: cd BuddySuite/buddysuite
    $: bash run_tests.sh

The second method should be run just before submitting a pull request, and uses
 [Docker](https://docs.docker.com/engine/installation/) to build a clean Linux environment with all dependencies. 
 The tests are then run in a container that mimics the environment used by Travis-CI to monitor the state of the BuddySuite
 repository:
    
    $: cd BuddySuite
    $: docker build -t docker-build:latest docker-build
    $: docker run -v $PWD:/home/docker/BuddySuite docker-build:latest
 
 Or just run the provided shell script:
 
    $: ./run_docker.sh

## Citation
There is a very short application note on bioRxiv that can be cited if you use BuddySuite in your work.

[DOI: 10.1101/040675](http://dx.doi.org/10.1101/040675)


## Contact
Any comments you may have would be really appreciated. Please feel free to add issues in the GitHub issue tracker or
 contact Steve Bond (lead developer) directly at [steve.bond@nih.gov](mailto:steve.bond@nih.gov)
