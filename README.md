# BuddySuite
Do fun stuff with biological data files. Seriously, biological data is fun stuff :)

## Description
The BuddySuite modules are designed to be 'one-stop-shop' command line tools for common biological data file 
manipulations.

Currently the Buddy 'Suite' only consists of one stable module, 
[SeqBuddy](https://github.com/biologyguy/BuddySuite/wiki/SeqBuddy), although three more modules are currently under 
development using the same general architecture:

- [AlignBuddy](https://github.com/biologyguy/BuddySuite/wiki/AlignBuddy)
- [DatabaseBuddy](https://github.com/biologyguy/BuddySuite/wiki/DatabaseBuddy)
- [PhyloBuddy](https://github.com/biologyguy/BuddySuite/wiki/PhyloBuddy)

Being pure Python, the BuddySuite should be cross platform. Almost all development and testing has been done on Linux
  and Mac OS X, however, so if you are a Windows user experiencing weird behavior, please let me know.

## Standalone installation
A graphical installer is in the works, but for now simply download the desired Buddy tool(s), and make executable.
    
    $: wget https://raw.github.com/biologyguy/BuddySuite/master/SeqBuddy
    $: chmod +x SeqBuddy

Run with the -h flag to see a list of available functions.

    $: ./SeqBuddy -h

I like to sym-link the main tools to short commands somewhere in my PATH (e.g., 'sb' for SeqBuddy, 'alb' for AlignBuddy, 
etc.), so they can be accessed quickly (examples in the wiki use these short forms).

    $: ln -s /path/to/SeqBuddy.py /usr/local/bin/sb
   
These are the stable release versions of the BuddySuite. If bugs are found they will be fixed, but the *expected* 
behavior will not be changed once the release is finalized. Likewise, new features added to the development versions
will not become available in the standalones until the next release. Version can be displayed using the -v flag.

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
[here](https://www.python.org/downloads/). 

The blast, bl2seq, and purge functions require access to the blastp, blastn, and blastdbcmd binaries from the 
[NCBI C++ toolkit](http://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/).
 
[BioPython](http://biopython.org/) is used heavily by the suite, and the package must be installed to use the development
version of the software. Furthermore, any BioPython versions earlier than 16.6 will cause unit tests to fail. 
This dependency is bundled into the stand alone version, however, so no extra download is required unless you are 
developing (or want the bleeding edge).

## [SeqBuddy](https://github.com/biologyguy/BuddySuite/wiki/SeqBuddy)
### Modifying flags
*Flag* | *Description*
------ | ----------
-a --alphabet | Sets the alphabet for the sequence. Otherwise, alphabet will be guessed.
-i --in_place | Rewrites the input file in-place. Be careful!
-o --out_format | Specify the format you want the output returned in
-p --params | Some functions can be uniquely modified by -p; see function for details
-q --quiet | Suppress stderr messages
-v --version | Output version information

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
[pull_random_record](https://github.com/biologyguy/BuddySuite/wiki/Pull-random-record) | -prr | None | Extract random sequences. Use the -p flag to increase the number of sequences returned
[pull_record_ends](https://github.com/biologyguy/BuddySuite/wiki/Rename-IDs) | -pre | \<amount (int)\> \<'front'\|'rear'\> | Get the ends (front or rear) of all sequences in a file
[extract_region](https://github.com/biologyguy/BuddySuite/wiki/Extract-region) | -er | \<start (int)\> \<end (int)\> | Pull out sub-sequences
[delete_records](https://github.com/biologyguy/BuddySuite/wiki/Delete-records) | -dr | \<regex pattern(s)\> | Remove reocrds from a file. The deleted IDs are sent to stderr
[delete_features](https://github.com/biologyguy/BuddySuite/wiki/Delete-features) | -df | \<regex pattern(s)\> | Remove specified features from all records
[delete_large](https://github.com/biologyguy/BuddySuite/wiki/Delete-large) | -dsm | \<threshold (int)\> | Delete sequences with length below threshold
[delete_small](https://github.com/biologyguy/BuddySuite/wiki/Delete-small) | -dlg | \<threshold (int)\> | Delete sequences with length above threshold
[find_repeats](https://github.com/biologyguy/BuddySuite/wiki/Find-repeats) | -frp | None | Identify whether a file contains repeat sequences and/or sequence ids
[delete_repeats](https://github.com/biologyguy/BuddySuite/wiki/Delete-repeats) | -drp | None | Strip repeat records (ids and/or identical sequences)
[merge](https://github.com/biologyguy/BuddySuite/wiki/Merge) | -mg | None | Group a bunch of seq files together
[blast](https://github.com/biologyguy/BuddySuite/wiki/BLAST) | -bl | \<BLAST database\> | BLAST your sequence file using common blast settings, return the hits from blastdb
[bl2seq](https://github.com/biologyguy/BuddySuite/wiki/Blast-2-seqs) | -bl2s | None | All-by-all blast among sequences using bl2seq. Only Returns top hit from each search
[purge](https://github.com/biologyguy/BuddySuite/wiki/Purge) | -prg | \<Max BLAST bit-score (int)\> | Delete sequences with high similarity
[guess_alphabet](https://github.com/biologyguy/BuddySuite/wiki/Guess-alphabet) | -ga | None | Return the alphabet type found in the input file
[guess_format](https://github.com/biologyguy/BuddySuite/wiki/Guess-format) | -gf | None | Guess the flatfile format of the input file
[molecular_weight](https://github.com/biologyguy/BuddySuite/wiki/Molecular-weight) | -mw | None | Computes the molecular weight of all of the sequences found in the input file.
[split_by_taxa](https://github.com/biologyguy/BuddySuite/wiki/Split-by-taxa) | -sbt | \<split char (s)\>, \<out dir (s)\> | Removes all sequences that do not contain a specific taxa in their ID.
[isoelectric_point](https://github.com/biologyguy/BuddySuite/wiki/Isoelectric-point) | -ip | None | Returns a list of isoelectric points for each protein sequence in the file.
[count_residues](https://github.com/biologyguy/BuddySuite/wiki/Count-residues) | -cr | None | Returns a list of residues and their frequencies for each sequence in the file.
[find_restriction_sites](https://github.com/biologyguy/BuddySuite/wiki/Find-restriction-sites) | -sbt | \<commercial (s)\>, \<num_cuts (int)\>| Returns a dictionary of all of the restriction sites and their indices for each sequence in the file. 
