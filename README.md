<p align="center"><a href="https://github.com/biologyguy/BuddySuite/wiki">
<img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/images/BuddySuite-logo.gif" /></a></p>
<p align="center">
<a href="https://github.com/biologyguy/BuddySuite/wiki/SeqBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/images/SeqBuddy-logo.gif" width=20%/></a>
<a href="https://github.com/biologyguy/BuddySuite/wiki/AlignBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/images/AlignBuddy-logo.gif" width=25%/></a>
<a href="https://github.com/biologyguy/BuddySuite/wiki/DBBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/images/DBBuddy-logo.gif" width=20%/></a>
<a href="https://github.com/biologyguy/BuddySuite/wiki/PhyloBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/images/PhyloBuddy-logo.gif" width=25%/></a>
</p>
<p align="center">Do fun stuff with biological data files. Seriously, biological data is fun stuff :)</p>
___
## Description
The BuddySuite modules are designed to be 'one-stop-shop' command line tools for common biological data file 
manipulations.

Currently the Buddy 'Suite' only consists of two relatively mature modules, 
[SeqBuddy](https://github.com/biologyguy/BuddySuite/wiki/SeqBuddy) and 
[AlignBuddy](https://github.com/biologyguy/BuddySuite/wiki/AlignBuddy), although two more modules are currently under 
development using the same general architecture 
([DatabaseBuddy](https://github.com/biologyguy/BuddySuite/wiki/DatabaseBuddy) and 
[PhyloBuddy](https://github.com/biologyguy/BuddySuite/wiki/PhyloBuddy))

Being pure Python, the BuddySuite should be cross platform. Almost all development and testing has been done on Linux
  and Mac OS X, however, so if you are a Windows user experiencing weird behavior, please let me know.

## Standalone installation 
### This is still an Alpha version, but it seems to be working for Mac and Linux
Download the graphical installer and run from the command line
    
    $: wget https://raw.github.com/biologyguy/BuddySuite/master/BuddySuite.py
    $: chmod +x BuddySuite.py
    $: ./BuddySuite.py

By default, the installer will place short form sym-links to the main tools in your PATH (e.g., 'sb' for SeqBuddy, 'alb'
 for AlignBuddy, etc.), so they can be accessed quickly (examples in the wiki use these short forms). The full names of
 each tool will also be added to PATH.

The installer bundles stable release versions of the BuddySuite. If bugs are found they will be fixed, but the *expected* 
behavior will not be changed once the release is finalized. Likewise, new features added to the development versions
will not become available in the installer until the next release. Versions of each tool or the installer can be 
displayed using the -v flag.

## Development version installation and contribution
All new features are developed in the 'workshop' versions of the buddy programs. These may be less stable than the 
official release versions, and have dependencies (see below).
The easiest way to get the development version up and running is to clone the repository.

    $: git clone https://github.com/biologyguy/BuddySuite.git

The Buddy tools are structured so they can be used as importable modules as well as command line programs. If you wish
to contribute to the project, new features require three components:

1. A self contained function that accepts a buddy object as input, and (usually) returns a new buddy object.
2. An argparse entry, allowing the function to be called from the command line.
3. Wrapper code in the `if __name__ == '__main__':` block to handle command line calls


## Dependencies
The BuddySuite is written in Python3 and is not backwards compatible with Python2. Python3 can be downloaded from 
[here](https://www.python.org/downloads/), or use (Anaconda)[http://continuum.io/downloads#py34] as a more comprehensive 
solution. 

The SeqBuddy blast, bl2seq, and purge functions require access to the blastp, blastn, and blastdbcmd binaries from the 
[NCBI C++ toolkit](http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/).
 
[BioPython](http://biopython.org/) is used heavily by the suite, and the package must be installed to use the development
version of the software. Furthermore, any BioPython versions earlier than 16.6 will cause unit tests to fail. 
This dependency is bundled into the stand alone version, however, so no extra download is required unless you are 
developing (or want the bleeding edge).

AlignBuddy and PhyloBuddy will eventually wrap a number of third party tools, but similar to the BLAST binaries, they
 will only be necessary to use their respective functions.