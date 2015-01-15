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
 if you want to use the blast, bl2seq, or purge functions
 

## SeqBuddy functions
*Function* | *Flag* | *Parameters* | *Description*
---------- | -------- | ---------- | ----------
clean_seq | -cs | None | Strip out non-sequence characters, such as stops (*) and gaps (-)
uppercase | -uc | None | Convert all sequences to uppercase
lowercase | -lc | None | Convert all sequences to lowercase
delete_metadata | -dm | None | Remove meta-data from file (only id is retained)
raw_seq | -rs | None | Return line break separated sequences
translate | -tr | None | Convert coding sequences into amino acid sequences
translate6frames | -tr6 | None | Translate nucleotide sequences into all six reading frames
back_translate | -btr | None | Convert amino acid sequences into codons. Select mode with -p flag ['random', <others>]
transcribe | -d2r | None | Convert DNA sequences to RNA
back_transcribe | -r2d | None | Convert RNA sequences to DNA
list_ids | -li | None | Output all the sequence identifiers in a file. Use -p to specify # columns to write
num_seqs | -ns | None | Counts how many sequences are present in an input file
concat_seqs | -cts | None | Concatenate a bunch of sequences into a single solid string
map_features_dna2prot | -fd2p | None | Take the features annotated onto cDNA sequences and map to protein sequences. Both a protein and cDNA file must be passed in
map_features_prot2dna | -fp2d | None | Take the features annotated onto protein sequences and map to cDNA sequences. Both a protein and cDNA file must be passed in
rename_ids | -ri | <regex pattern> <subst string> | Replace some pattern in ids with something else.
combine_features | -cf | None | Takes the features in two files and combines them for each sequence
order_features_by_position | -ofp | None | Change the output order of sequence features, based on sequence position
order_features_alphabetically | -ofa | None | Change the output order of sequence features, based on sequence position
screw_formats | -sf | <new format> | Change the file format to something else
shuffle | -sh | None | Randomly reorder the position of records in the file
hash_seq_ids | -hsi | None | Rename all the identifiers in a sequence list to a 10 character hash
pull_records | -pr | <regex pattern> | Get all the records with ids containing a given string
pull_record_ends | -pre | <amount (int)> <'front'|'rear'> | Get the ends (front or rear) of all sequences in a file
delete_records | -dr | <regex pattern(s)> | Remove reocrds from a file. The deleted IDs are sent to stderr
delete_features | -df | <regex pattern(s)> | Remove specified features from all records
delete_repeats | -drp | None | Strip repeat records (ids and/or identical sequences
find_repeats | -fr | None | Identify whether a file contains repeat sequences and/or sequence ids
merge | -mg | None | Group a bunch of seq files together
blast | -bl | <BLAST database> | BLAST your sequence file using common blast settings, return the hits from blastdb
bl2seq | -bl2s | None | All-by-all blast among sequences using bl2seq. Only Returns top hit from each search
purge | -prg | <Max BLAST score (int)> | Delete sequences with high similarity
guess_alphabet | -ga | None | Return the alphabet type found in the input file
guess_format | -gf | None | Guess the flatfile format of the input file