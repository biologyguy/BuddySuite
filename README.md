##This is a branch of the main project.
Please head to the [main](https://github.com/biologyguy/BuddySuite) page for the latest working version.

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

Being pure Python, the BuddySuite should be cross platform. Development and testing have been done on Linux
 and Mac OS X, however, so it is unclear if the Suite will work correctly under Windows.

## Dependencies
The BuddySuite is written in Python3 and is not backwards compatible with Python2. Python3 can be downloaded from
 [here](https://www.python.org/downloads/), or use [Anaconda](http://continuum.io/downloads#py34) as a more
 comprehensive solution. 

The SeqBuddy blast, bl2seq, and purge functions require access to the blastp, blastn, and blastdbcmd binaries from the
 [NCBI C++ toolkit](http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/). If not already in your PATH the binaries will be
 downloaded by the installer. SeqBuddy.py will also attempt to download the binaries if any BLAST dependant functions
 are called and the programs are not found in PATH.
 
[BioPython](http://biopython.org/) is used heavily by the suite, and the package must be installed to use the
 development version of the software. Furthermore, any BioPython versions earlier than 16.6 will cause unit tests to
 fail. This dependency is bundled with the installer, however, so no extra download is required unless you are
 developing (or want the bleeding edge).

AlignBuddy and PhyloBuddy will eventually wrap a number of third party tools, but similar to the BLAST binaries, they
 will only be necessary to use their respective functions.
 
## Standalone installation 
#### This is still an Alpha version, but it seems to be working for Mac and Linux. It will not work on Windows.
Download the graphical installer and run it from the command line
    
    $: wget https://raw.github.com/biologyguy/BuddySuite/master/BuddySuite.py --no-check-certificate
    $: chmod +x BuddySuite.py
    $: ./BuddySuite.py

By default, the installer will place short form sym-links to the main tools in your PATH (e.g., 'sb' for SeqBuddy, 'alb'
 for AlignBuddy, etc.), so they can be accessed quickly (examples in the wiki use these short forms). The full names of
 each tool will also be added to PATH. If working outside the context of a graphical OS (on a cluster, for example), the
 installer may be run with the -cmd flag, which will walk you through the install processes directly from the command
 line.

Once out of alpha, the installer will only bundle stable release versions of the BuddySuite. If bugs are found they will
 be fixed, but the *expected* behavior will not be changed once the release is finalized. Likewise, new features added
 to the development versions will not become available in the installer until the next release. Versions of each tool or
 the installer can be displayed using the -v flag.

## Development version installation and contribution
All new features are developed in the 'workshop' versions of the buddy programs. These may be less stable than the
 official release versions, and may have extra dependencies. The easiest way to get the development version
 up and running is to [clone/fork](https://help.github.com/articles/fork-a-repo/) the repository.

    $: git clone https://github.com/biologyguy/BuddySuite.git

The Buddy tools are structured so they can be used as importable modules as well as command line programs. If you wish
to contribute to the project, new features require three components:

1. A self contained function that accepts a buddy object as input, and (usually) returns a new buddy object.
2. An argparse entry, allowing the function to be called from the command line.
3. Wrapper code in the `if __name__ == '__main__':` block to handle command line calls

The Buddy classes (i.e., SeqBuddy, AlignBuddy, etc) contain a number of core attributes and methods (see below), and act
 as a standardized input/output for all functions in the Suite. Instantiation of a buddy object is designed to be
 extremely flexible, as the detection and handling of most input types (stdin, file handles, file paths, list of
 SeqRecords, or plain text) are all cooked into the classes.

##### SeqBuddy
    sb_obj = SeqBuddy(_input, _in_format=None, _out_format=None, _alpha=None)

*Attribute* | *Description*
----------- | -------------
alpha | An IUPAC object from Bio.Alphabet (one of: IUPAC.protein, IUPAC.ambiguous_dna, or IUPAC.ambiguous_rna). Plain text representatives of each (e.g., 'dna', 'prot', 'protein', 'r') will be understood by the SeqBuddy \__init\__() method, or the alphabet will be guessed if not explicitly set.
in_format | The [flat file format](http://biopython.org/wiki/SeqIO#File_Formats) sequences are read from. If explicitly set, SeqBuddy will only attempt to read the file in the given format (returning no sequences if the wrong format is specified), otherwise it will guess the format.
out_format | Controls the format used when SeqBuddy objects are written. By default, this will be the same as in_format.
records | A list of Bio.SeqRecord objects.
 
*Method* | *Description*
-------- | -------------
to_dict() | Return a dictionary of everything in self.records using SeqRecord.id as keys
print() | Write all records to stdout using out_format
write(_file_path) | Write all records to file using out_format

##### AlignBuddy
    alb_obj = AlignBuddy(_input, _in_format=None, _out_format=None)

*Attribute* | *Description*
----------- | -------------
alpha | An IUPAC object from Bio.Alphabet, same as in SeqBuddy. The constructor does not accept this explicitly, it is guessed from the sequences in the alignment(s).
in_format | The [flat file format](http://biopython.org/wiki/AlignIO#File_Formats) sequences are read from. If explicitly set, AlignBuddy will only attempt to read the file in the given format (returning no alignments if the wrong format is specified), otherwise it will guess the format.
out_format | Controls the format used when AlignBuddy objects are written. By default, this will be the same as in_format.
alignments | A list of Bio.Align objects.
 
*Method* | *Description*
-------- | -------------
print() | Write all alignments to stdout using out_format
write(_file_path) | Write all alignments to file using out_format

##### DbBuddy
    dbb_obj = DbBuddy(_input, _database=None, _out_format="summary")

Note: The core of DbBuddy is still under active development, so expect weirdness.

*Attribute* | *Description*
----------- | -------------
search_terms | List of search terms used to query public databases
records | Dictionary of Record objects (not BioPython records!), using accession numbers as keys.
recycle_bin | Also a dictionary of records, but which have been filtered out of the main records dict (i.e., will not be output by print())
out_format | Controls the format used when DbBuddy objects are written. Valid formats include "summary", "full_summary", "ids", "accessions", and any of the supported [BioPython SeqIO](http://biopython.org/wiki/SeqIO#File_Formats) formats.
failures | A list of Failure objects. The Failure class is used to track issues encountered while communicating with the public databases.
databases | A list of databases that DbBuddy will query. Valid options include "all", "uniprot", "ensemble", "gb", and "refseq".
 
*Method* | *Description*
-------- | -------------
record_breakdown() | Return a dictionary with counts for 'full', 'partial', and 'accession' only records.
filter_records(regex) | Send any records not matched by 'regex' to 'recycle_bin'
restore_records(regex) | Send any records matching 'regex' from 'recycle_bin' back to 'records'
print(_num=0, quiet=False, columns=None, destination=None) | Write records to stdout or a path (set with \_destination). The number of records and the columns displayed can be set with their respective arguments.

##### PhyloBuddy
    pb_obj = PhyloBuddy(_input)

PhyloBuddy is not implemented yet.
 
 
## Contact
Any comments you may have would be really appreciated. Please feel free to add issues in the GitHub issue tracker or
 contact me directly at [steve.bond@nih.gov](mailto:steve.bond@nih.gov)