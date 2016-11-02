#!/usr/bin/env python3

import sys
import argparse
import timeit

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="performanceScanner", description="Check function time", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("regex", help="pattern to search", action="store")
    # parser.add_argument("replace", help="replacement", action="store")
    # parser.add_argument("file", help="input file", action="store")
    parser.add_argument("-m", "--module", help="specify a single module", action='store', default = "seqbuddy")
    parser.add_argument("-c", "--command", help="specify a single module", action='store')
    parser.add_argument("-i", "--iteration", help="specify iteration number", action='store', default=10)
    # use seq buddy for now, then ad the --module option
    in_args = parser.parse_args()

    # create a dictionary with the default parameters for each function of seqbuddy (-c, --command)
    opts_sb = {
        "annotate": "misc_feature 1-10", "ave_seq_length": "", "back_translate": "", "bl2seq": "", "blast": "mnemiopsis_leidyi",
        "clean_seq": "", "complement": "", "concat_seqs": "", "count_codons": "", "count_residues": "concatenate",
        "degenerate_sequence": "1", "delete_features": "TMD1", "delete_large": "400", "delete_metadata": "",
        "delete_records": "Aae-", "delete_repeats": "", "delete_small": "100", "extract_feature_sequences": "TMD1",
        "extract_regions": ":50", "find_CpG": "-o genbank", "find_orfs": "-o genbank", "find_pattern": "\'LLL\' \'LLY\'",
        "find_repeats": "3", "find_restriction_sites": "commercial", "group_by_prefix": "~/BuddySuite/diagnostics/files_prefix/",
        "group_by_regex": "~/BuddySuite/diagnostics/files_regex/ Sra", "guess_alphabet": "~/BuddySuite/diagnostics/files_regex/*",
        "guess_format": "~/BuddySuite/diagnostics/files_regex/*", "hash_seq_ids": "10", "insert_seq": "DYKDDDDK 0",
        "isoelectric_point": "", "list_features": "", "list_ids": "3", "lowercase": "", "make_ids_unique": "10 \'-\'",
        "map_features_nucl2prot": "All_pannexins_nuc.gb All_pannexins_pep.fa",
        "map_features_prot2nucl": "All_pannexins_nuc.fa All_pannexins_pep.gb", "merge": "All_pannexins_pep.gb All_pannexins_pep2.gb",
        "molecular_weight": "", "num_seqs": "", "order_features_alphabetically": "rev", "order_features_by_position": "rev",
        "order_ids": "rev", "order_ids_randomly": "", "prosite_scan": "strict", "pull_random_record": "4", "pull_records": "Aae-Panx",
        "pull_record_ends": "-100", "pull_records_with_feature": "PHOSPHO", "purge": "230", "rename_ids": "Mle Mnemiopsis",
        "replace_subseq": "IL].{1,4}[IL] _motif_", "reverse_complement": "", "reverse_transcribe": "", "screw_formats": "gb",
        "select_frame": "3", "shuffle_seqs": "", "translate": "", "translate6frames": "", "transcribe": "",
        "transmembrane_domains": "", "uppercase": "",
    }

    # create a dictionary with the default parameters for each function of alignbuddy (-c, --command)
    opts_ab = {
        "alignment_lengths": "", "back_transcribe": "", "bootstrap": "3", "clean_seq": "strict", "concat_alignments": "[a-z]{3}-Panx",
        "consensus": "", "delete_records": "PanxγA", "enforce_triplets": "", "extract_regions": "10 110", "generate_alignment": "mafft",
        "hash_ids": "10", "list_ids": "3", "lowercase": "", "mapfeat2align": "All_pannexins_pep.gb", "num_seqs": "",
        "order_ids": "rev", "pull_records": "PanxβB PanxβC", "rename_ids": "Mle Mnemiopsis", "screw_formats": "fasta",
        "split_to_files": "~/BuddySuite/diagnostics/files_al Pannexin_", "transcribe": "", "translate": "", "trimal": "clean", "uppercase": "",
    }

    # create a dictionary with the default parameters for each function of phylobuddy (-c, --command)
    # All_pannexins_pep_aln_tree.nex --> default file
    opts_pb = {
        "collapse_polytomies": "length 0.1", "consensus_tree": "0.5", "display_trees": "", "distance": "uwrf", "generate_tree": "raxmlHPC-SSE3 -o nexus",
        "hash_ids": "3", "list_ids": "3", "num_tips": "", "print_trees": "", "prune_taxa": "Cte-PanxζI", "rename_ids": "α A",
        "root": "2", "screw_formats": "newick", "show_unique": "", "split_polytomies": "", "unroot": "",
    }

    if in_args.module == "seqbuddy":
        if in_args.command == "all":
            for flag, args in opts_sb.items():
                sys.stdout.write("%s: " % flag)
                sys.stdout.flush()
                if flag in ["concat_seqs", "delete_metadata", "pull_records_with_feature", "delete_features", "extract_feature_sequences", "extract_regions"]:
                    timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.gb --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                    sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
                elif flag in ["degenerate_sequence", "find_CpG", "find_orfs", "find_restriction_sites", "translate", "translate6frames", "transcribe","uppercase"]:
                    timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_nuc.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                    sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
                elif flag in ["guess_alphabet", "guess_format", "map_features_nucl2prot", "map_features_prot2nucl", "merge"]:
                    timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy \'%s\' --%s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (args, flag), number=int(in_args.iteration))
                    sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
                elif flag == "reverse_transcribe":
                    timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_rna.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                    sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
                else:
                    timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                    sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
        elif in_args.command != "all" and in_args.command in opts_sb:
            sys.stdout.write("%s: " % in_args.command)
            sys.stdout.flush()
            if in_args.command in ["concat_seqs", "delete_metadata", "pull_records_with_feature", "delete_features", "extract_feature_sequences", "extract_regions"]:
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.gb --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_sb[in_args.command]), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
            elif in_args.command in ["degenerate_sequence", "find_CpG", "find_orfs", "find_restriction_sites", "translate", "translate6frames", "transcribe", "uppercase"]:
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_nuc.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_sb[in_args.command]), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
            elif in_args.command in ["guess_alphabet", "guess_format", "map_features_nucl2prot", "map_features_prot2nucl", "merge"]:
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy \'%s\' --%s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (opts_sb[in_args.command], in_args.command), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
            elif in_args.command == "reverse_transcribe":
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_rna.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_sb[in_args.command]), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
            else:
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_sb[in_args.command]), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
        elif in_args.command not in opts_sb:
            sys.stdout.write("Invalid command\n")
    elif in_args.module == "alignbuddy":
        if in_args.command == "all":
            for flag, args in opts_ab.items():
                sys.stdout.write("%s: " % flag)
                sys.stdout.flush()
                if flag == "generate_alignment":
                    timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("alignbuddy All_pannexins_pep.gb --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                    sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
                elif flag == "transcribe":
                    timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("alignbuddy All_pannexins_nuc_aln.gb --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                    sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
                elif flag == "translate":
                    timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("alignbuddy All_pannexins_rna_aln.gb --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                    sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
                else:
                    timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("alignbuddy All_pannexins_pep_aln.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                    sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
        elif in_args.command != "all" and in_args.command in opts_ab:
            sys.stdout.write("%s: " % in_args.command)
            sys.stdout.flush()
            if in_args.command == "generate_alignment":
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("alignbuddy All_pannexins_pep.gb --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_ab[in_args.command]), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
            elif in_args.command == "transcribe":
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("alignbuddy All_pannexins_nuc_aln.gb --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_ab[in_args.command]), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
            elif in_args.command == "translate":
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("alignbuddy All_pannexins_rna_aln.gb --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_ab[in_args.command]), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
            else:
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("alignbuddy All_pannexins_pep_aln.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_ab[in_args.command]), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
        elif in_args.command not in opts_ab:
            sys.stdout.write("Invalid command\n")
    elif in_args.module == "phylobuddy":
        if in_args.command == "all":
            for flag, args in opts_pb.items():
                sys.stdout.write("%s: " % flag)
                sys.stdout.flush()
                if flag == "consensus_tree":
                    timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("phylobuddy All_pannexins_pep.nwk --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                    sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
                elif flag == "generate_tree":
                    timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("phylobuddy All_pannexins_pep_aln.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                    sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
                else:
                    timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("phylobuddy All_pannexins_pep_aln_tree.nex --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                    sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
        elif in_args.command != "all" and in_args.command in opts_pb:
            sys.stdout.write("%s: " % in_args.command)
            sys.stdout.flush()
            if in_args.command == "consensus_tree":
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("phylobuddy All_pannexins_pep.nwk --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_pb[in_args.command]), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
            elif in_args.command == "generate_tree":
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("phylobuddy All_pannexins_pep_aln.gb --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_pb[in_args.command]), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
            else:
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("phylobuddy All_pannexins_pep_aln_tree.nex --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_pb[in_args.command]), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
        else:
            sys.stdout.write("Invalid command\n")
    else:
        sys.stdout.write("Invalid module\n")