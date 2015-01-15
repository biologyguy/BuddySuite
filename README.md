# BuddySuite
Do fun stuff with biological data files. Seriously, biological data is fun stuff :)

## Description
The BuddySuite modules are designed to be 'one-stop-shop' command line tools for common biological data file 
manipulations.

Currently the Buddy'Suite' only consists of one module, SeqBuddy, with plans to create at least three more 
modules using the same general architecture:

- AlignBuddy
- DatabaseBuddy
- PhyloBuddy.

## Installation
Simply download SeqBuddy.py and make it executable
 
 `chmod +x SeqBuddy.py`

Run with the -h flag to see a list of available functions

  `./SeqBuddy.py -h`
  
## Dependencies
SeqBuddy requires the [BioPython](http://biopython.org/) package

You will need blastp, blastn, and blastdbcmd from the [NCBI C++ toolkit](http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/)
 if you want to use the blast or bl2seq functions
 

## SeqBuddy functions
*function* | *Flag* | *Parameters* | *Description*
---------- | ------ | ---------- | ----------
clean_seq | cs | None | Strip out non-sequence characters, such as stops () and gaps ()
 