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

I like to sym-link SeqBuddy to the command 'sb' somewhere in my PATH, so I can access it quickly. For example:
 
 `ln -s /path/to/SeqBuddy.py /usr/local/bin/sb`


## Dependencies
SeqBuddy is written in Python3 and is not backwards compatible with Python2. Python3 can be downloaded from 
[here](https://www.python.org/downloads/) 

The [BioPython](http://biopython.org/) package is also required. 

You will need blastp, blastn, and blastdbcmd from the [NCBI C++ toolkit](http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/)
 if you want to use the blast, bl2seq, or purge functions
 

## [SeqBuddy](https://github.com/biologyguy/BuddySuite/wiki/SeqBuddy)
### Modifying flags
*Flag* | *Description*
------ | ----------
-o --out_format | Specify the format you want the output returned in
-i --in_place | Rewrites the input file in-place. Be careful!
-p --params | Some functions can be uniquely modified by -p; see function for details
-q --quiet | Suppress stderr messages

### Functions
*Function* | *Flag* | *Parameters* | *Description*
---------- | -------- | ---------- | ----------
[clean_seq](https://github.com/biologyguy/BuddySuite/wiki/Clean-sequence) | -cs | None | Strip out non-sequence characters, such as stops (*) and gaps (-)
[uppercase](https://github.com/biologyguy/BuddySuite/wiki/Uppercase) | -uc | None | Convert all sequences to uppercase
[lowercase](https://github.com/biologyguy/BuddySuite/wiki/Lowercase) | -lc | None | Convert all sequences to lowercase
[delete_metadata](https://github.com/biologyguy/BuddySuite/wiki/Delete-metadata) | -dm | None | Remove meta-data from file (only id is retained)
[raw_seq](https://github.com/biologyguy/BuddySuite/wiki/Raw-sequence) | -rs | None | Return line break separated sequences
[translate](https://github.com/biologyguy/BuddySuite/wiki/Translate) | -tr | None | Convert coding sequences into amino acid sequences
[select_frame](https://github.com/biologyguy/BuddySuite/wiki/Select-frame) | -sfr | \<frame (1, 2, or 3)\> | Change the reading from of sequences by deleting characters off of the front
[translate6frames](https://github.com/biologyguy/BuddySuite/wiki/Translate-6-frames) | -tr6 | None | Translate nucleotide sequences into all six reading frames
[back_translate](https://github.com/biologyguy/BuddySuite/wiki/Back-translate) | -btr | None | Convert amino acid sequences into codons. Select mode/species with -p flag \['random', 'optimized'\] \['human', 'mouse', 'yeast', 'ecoli'\]
[transcribe](https://github.com/biologyguy/BuddySuite/wiki/Transcribe) | -d2r | None | Convert DNA sequences to RNA
[back_transcribe](https://github.com/biologyguy/BuddySuite/wiki/Back-transcribe) | -r2d | None | Convert RNA sequences to DNA
[complement](https://github.com/biologyguy/BuddySuite/wiki/Complement) | -cmp | None | Return complement of nucleotide sequence
[reverse_complement](https://github.com/biologyguy/BuddySuite/wiki/Reverse-complement) | -rc | None | Return reverse complement of nucleotide sequence
[list_ids](https://github.com/biologyguy/BuddySuite/wiki/List-IDs) | -li | None | Output all the sequence identifiers in a file. Use -p to specify # columns to write
[num_seqs](https://github.com/biologyguy/BuddySuite/wiki/Number-of-sequences) | -ns | None | Counts how many sequences are present in an input file
[ave_seq_length](https://github.com/biologyguy/BuddySuite/wiki/Average-sequence-length) | -asl | None | Find the average length of all sequences in an input file
[concat_seqs](https://github.com/biologyguy/BuddySuite/wiki/Concatinate-sequences) | -cts | None | Concatenate a bunch of sequences into a single solid string
[map_features_dna2prot](https://github.com/biologyguy/BuddySuite/wiki/Map-features-nucl2prot) | -fd2p | None | Take the features annotated onto cDNA sequences and map to protein sequences. Both a protein and cDNA file must be passed in
[map_features_prot2dna](https://github.com/biologyguy/BuddySuite/wiki/Map-features-prot2nucl) | -fp2d | None | Take the features annotated onto protein sequences and map to cDNA sequences. Both a protein and cDNA file must be passed in
[rename_ids](https://github.com/biologyguy/BuddySuite/wiki/Rename-IDs) | -ri | \<regex pattern\> \<subst string\> | Replace some pattern in ids with something else.
[combine_features](https://github.com/biologyguy/BuddySuite/wiki/Combine-features) | -cf | None | Takes the features in two files and combines them for each sequence
[order_features_by_position](https://github.com/biologyguy/BuddySuite/wiki/Order-features-by-position) | -ofp | None | Change the output order of sequence features, based on sequence position
[order_features_alphabetically](https://github.com/biologyguy/BuddySuite/wiki/Order-features-alphabetically) | -ofa | None | Change the output order of sequence features, based on sequence position
[screw_formats](https://github.com/biologyguy/BuddySuite/wiki/Screw-formats) | -sf | \<new format\> | Change the file format to something else
[shuffle](https://github.com/biologyguy/BuddySuite/wiki/Shuffle-sequences) | -sh | None | Randomly reorder the position of records in the file
[order_ids](https://github.com/biologyguy/BuddySuite/wiki/Order-IDs) | -oi | None | Sort all sequences by id in alpha-numeric order. Use -p 'rev' for reverse order.
[hash_seq_ids](https://github.com/biologyguy/BuddySuite/wiki/Hash-sequence-IDs) | -hsi | None | Rename all the identifiers in a sequence list to a 10 character hash
[pull_records](https://github.com/biologyguy/BuddySuite/wiki/Pull-records) | -pr | \<regex pattern\> | Get all the records with ids containing a given string
[pull_record_ends](https://github.com/biologyguy/BuddySuite/wiki/Rename-IDs) | -pre | \<amount (int)\> \<'front'\|'rear'\> | Get the ends (front or rear) of all sequences in a file
[extract_region](https://github.com/biologyguy/BuddySuite/wiki/Extract-region) | -er | \<start (int)\> \<end (int)\> | Pull out sub-sequences
[delete_records](https://github.com/biologyguy/BuddySuite/wiki/Delete-records) | -dr | \<regex pattern(s)\> | Remove reocrds from a file. The deleted IDs are sent to stderr
[delete_features](https://github.com/biologyguy/BuddySuite/wiki/Delete-features) | -df | \<regex pattern(s)\> | Remove specified features from all records
[find_repeats](https://github.com/biologyguy/BuddySuite/wiki/Find-repeats) | -frp | None | Identify whether a file contains repeat sequences and/or sequence ids
[delete_repeats](https://github.com/biologyguy/BuddySuite/wiki/Delete-repeats) | -drp | None | Strip repeat records (ids and/or identical sequences)
[merge](https://github.com/biologyguy/BuddySuite/wiki/Merge) | -mg | None | Group a bunch of seq files together
[blast](https://github.com/biologyguy/BuddySuite/wiki/BLAST) | -bl | \<BLAST database\> | BLAST your sequence file using common blast settings, return the hits from blastdb
bl2seq | -bl2s | None | All-by-all blast among sequences using bl2seq. Only Returns top hit from each search
purge | -prg | \<Max BLAST score (int)\> | Delete sequences with high similarity
guess_alphabet | -ga | None | Return the alphabet type found in the input file
guess_format | -gf | None | Guess the flatfile format of the input file