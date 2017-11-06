[![Build Status](https://travis-ci.org/biologyguy/BuddySuite.svg?branch=master)](https://travis-ci.org/biologyguy/BuddySuite)
[![Coverage Status](https://img.shields.io/coveralls/biologyguy/BuddySuite/master.svg)](https://coveralls.io/github/biologyguy/BuddySuite?branch=master)
[![PyPi version](https://img.shields.io/pypi/v/buddysuite.svg)](https://pypi.python.org/pypi/buddysuite)
<p align="center"><a href="https://github.com/biologyguy/BuddySuite/wiki">
<img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/buddysuite/images/BuddySuite-logo.png" width=70%/></a></p>
<p align="center">
<a href="https://github.com/biologyguy/BuddySuite/wiki/SeqBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/buddysuite/images/SeqBuddy-logo.gif" width=20%/></a>
<a href="https://github.com/biologyguy/BuddySuite/wiki/AlignBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/buddysuite/images/AlignBuddy-logo.gif" width=25%/></a>
<a href="https://github.com/biologyguy/BuddySuite/wiki/DatabaseBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/buddysuite/images/DBBuddy-logo.gif" width=20%/></a>
<a href="https://github.com/biologyguy/BuddySuite/wiki/PhyloBuddy"><img src="https://raw.githubusercontent.com/biologyguy/BuddySuite/master/buddysuite/images/PhyloBuddy-logo.gif" width=25%/></a>
</p>
<p align="center">Do fun stuff with biological data files. Seriously, biological data is fun stuff :)</p>

----

## Description
The BuddySuite modules are 'one-stop-shop' command-line tools for common biological data file
 manipulations. Formats are detected automatically, conversions are seamless, and you can pipe into
 or out of the modules to build custom bioinformatics workflows, allowing you to spend more time analyzing 
 your sequences, alignments, and phylogenetic trees, instead of wrangling them.

For example, the following command reads in three sequence files (all in different formats), pulls out records with RefSeq identifiers,
 calls MAFFT to generate an alignment, shifts gaps to force a codon alignment, calls RAxML to infer a phylogeny, and then roots
 the tree at its midpoint.
 
`$: ﻿seqbuddy seqs1.gb seqs2.embl seqs3.fasta --pull_records "[XN]M" | alignbuddy --generate_alignment mafft | alignbuddy --enforce_triplets | phylobuddy --generate_tree raxmlHPC-SSE3 | phylobuddy --root`

BuddySuite is a Python3 project, developed and extensively tested on Linux and Mac OS X. Full release versions
 are also tested on Windows before release, so should work as expected on Vista and above.

## Getting started
The simplest way to get up and running is:

```bash
$: pip install buddysuite 
$: buddysuite -setup
```

Further instructions are available in the [installation guide](https://github.com/biologyguy/BuddySuite/wiki/Installation-Guide).

There is also a [Beginners' Guide](https://github.com/biologyguy/BuddySuite/wiki/Beginners-Guide) to show you the basics,
 as well as a more comprehensive [Tutorial](https://github.com/biologyguy/BuddySuite/wiki/Tutorial).

Each tool in the BuddySuite has been extensively documented in the [wiki](https://github.com/biologyguy/BuddySuite/wiki),
 complete with worked examples and explanations for all arguments/options.
 
* [SeqBuddy](https://github.com/biologyguy/BuddySuite/wiki/SeqBuddy) 
* [AlignBuddy](https://github.com/biologyguy/BuddySuite/wiki/AlignBuddy)
* [PhyloBuddy](https://github.com/biologyguy/BuddySuite/wiki/PhyloBuddy)
* [DatabaseBuddy](https://github.com/biologyguy/BuddySuite/wiki/DatabaseBuddy)

## Developers
All of the individual Buddy toolkits are located in the 'buddysuite' directory and the 
 ['develop' branch](https://github.com/biologyguy/BuddySuite/tree/develop) is where all new features have been
 implemented. If you're interested in contributing, please refer to the
 [developer page](https://github.com/biologyguy/BuddySuite/wiki/Developers) for further information.

## Citation
[![doi](https://img.shields.io/badge/doi-10.1093/molbev/msx089-blue.svg?style=flat)](https://doi.org/10.1093/molbev/msx089)

If you use BuddySuite in your work, we would greatly appreciate references to our manuscript:

```bibtex
@article{doi:10.1093/molbev/msx089,
author = {Bond, Stephen R. and Keat, Karl E. and Barreira, Sofia N. and Baxevanis, Andreas D.},
title = {BuddySuite: Command-Line Toolkits for Manipulating Sequences, Alignments, and Phylogenetic Trees},
journal = {Molecular Biology and Evolution},
volume = {34},
number = {6},
pages = {1543-1546},
year = {2017},
doi = {10.1093/molbev/msx089},
URL = { + http://dx.doi.org/10.1093/molbev/msx089},
eprint = {/oup/backfile/content_public/journal/mbe/34/6/10.1093_molbev_msx089/1/msx089.pdf}
}
```

Through a series of unfortunate events and a surprisingly unresponsive publishing company, the manuscript is behind a
 paywall... Thank goodness for [pre-prints on bioRxiv](https://doi.org/10.1101/040675)!



## Contact
Any comments you have would be really appreciated. Please feel free to add issues in the GitHub issue tracker or
 contact Steve Bond (lead developer) directly at [steve.bond@nih.gov](mailto:steve.bond@nih.gov).
