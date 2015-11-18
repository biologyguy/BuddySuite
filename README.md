<p align="center"><a href="https://github.com/biologyguy/BuddySuite/wiki">
<img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/images/BuddySuite-logo.gif" /></a></p>
<p align="center">
<a href="https://github.com/biologyguy/BuddySuite/wiki/SeqBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/images/SeqBuddy-logo.gif" width=20%/></a>
<a href="https://github.com/biologyguy/BuddySuite/wiki/AlignBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/images/AlignBuddy-logo.gif" width=25%/></a>
<a href="https://github.com/biologyguy/BuddySuite/wiki/DatabaseBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/images/DBBuddy-logo.gif" width=20%/></a>
<a href="https://github.com/biologyguy/BuddySuite/wiki/PhyloBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/images/PhyloBuddy-logo.gif" width=25%/></a>
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
 [NCBI C++ toolkit](http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/). If not already in your PATH the binaries will be
 downloaded by the installer. SeqBuddy.py will also attempt to download the binaries if any BLAST dependant functions
 are called and the programs are not found in PATH. [BioPython](http://biopython.org/) is used heavily by the entire 
 suite; any version earlier than 16.6 will cause unit tests to fail. PhyloBuddy requires 
 [DendroPy](https://pythonhosted.org/DendroPy/) and version 3.0 (beta) of the
 [ETE toolkit](http://etetoolkit.org/download/).
 
## Standalone installation 
The installer will only run on Mac and Linux. If you would like to try the BuddySuite on Windows,
 you will need to install the development version (see below).

[Download the graphical installer](https://raw.github.com/biologyguy/BuddySuite/master/BuddySuite.py)
 and run it from the command line
    
    $: cd /path/to/download/folder
    $: chmod +x BuddySuite.py
    $: ./BuddySuite.py

By default, the installer will create short-form symbolic links for the main tools in your PATH ('sb' for SeqBuddy, 'alb'
 for AlignBuddy, 'pb' for PhyloBuddy, and 'db' for DatabaseBuddy), so they can be accessed quickly ([examples in the
 wiki](https://github.com/biologyguy/buddysuite/wiki) use these short forms). The full names of each tool will also be
 added to PATH. If working outside the context of a graphical OS (on a cluster, for example), the installer will run
 in command-line mode (also accessible with the -cmd flag on graphical systems, if you prefer that).

Once the BuddySuite moves out of beta, the installer will only bundle stable release versions of the BuddySuite. 
 If bugs are found they will be fixed, but the *expected* behavior will not be changed once the release is finalized. 
 Likewise, new features added to the development versions will not become available in the installer until the 
 next release. Versions of each tool or the installer can be displayed using the -v flag.

## Development version installation
All new features are developed in the 'workshop' versions of the buddy programs. These may be less stable than the
 official release versions, and may have extra dependencies. The easiest way to get the development version
 up and running is to [clone/fork](https://help.github.com/articles/fork-a-repo/) the repository.

    $: git clone https://github.com/biologyguy/BuddySuite.git

For further information on dependencies and how to contribute to the project, please see the
 [developer page](https://github.com/biologyguy/BuddySuite/wiki/Developers).

## Contact
Any comments you may have would be really appreciated. Please feel free to add issues in the GitHub issue tracker or
 contact me directly at [steve.bond@nih.gov](mailto:steve.bond@nih.gov)