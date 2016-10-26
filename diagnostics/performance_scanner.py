#!/usr/bin/env python3

import sys
import argparse
import timeit

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="performanceScanner", description="Check function time", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser.add_argument("regex", help="pattern to search", action="store")
    # parser.add_argument("replace", help="replacement", action="store")
    # parser.add_argument("file", help="input file", action="store")
    # parser.add_argument("-m", "--module", help="specify a single module", action='store')
    parser.add_argument("-c", "--command", help="specify a single module", action='store')
    parser.add_argument("-i", "--iteration", help="specify iteration number", action='store', default=10)
    # use seq buddy for now, then ad the --module option
    in_args = parser.parse_args()

    # create a dictionary with the default parameters for each function (-c, --command)
    opts_sb = {
        "ave_seq_length": "", "back_translate": "",
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
        "order_ids": "rev", "order_ids_randomly": "", "pull_random_record": "4", "pull_records": "Aae-Panx",
        "pull_record_ends": "-100", "pull_records_with_feature": "PHOSPHO", "rename_ids": "Mle Mnemiopsis",
        "replace_subseq": "IL].{1,4}[IL] _motif_",
    }

    # out of the dictionary for now (they take too long): "bl2seq": "", "blast": "mnemiopsis_leidyi", "prosite_scan": "strict", "purge": "230",

    if in_args.command == "all":
        for flag, args in opts_sb.items():
            sys.stdout.write("%s: " % flag)
            sys.stdout.flush()
            if flag == "concat_seqs" or flag == "delete_metadata" or flag == "pull_records_with_feature":
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.gb --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
            elif flag == "degenerate_sequence" or flag == "find_CpG" or flag == "find_orfs" or flag == "find_restriction_sites":
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_nuc.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
            elif flag == "delete_features" or flag == "extract_feature_sequences" or flag == "extract_regions":
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.gb --%s \'%s\'", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
            elif flag == "delete_records":
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.fa --%s \'%s\'", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
            elif flag == "guess_alphabet" or flag == "guess_format" or flag == "map_features_nucl2prot" or flag == "map_features_prot2nucl" or flag == "merge":
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy \'%s\' --%s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (args, flag), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
            else:
                timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (flag, args), number=int(in_args.iteration))
                sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
    elif in_args.command != "all" and in_args.command in opts_sb:
        sys.stdout.write("%s: " % in_args.command)
        sys.stdout.flush()
        if in_args.command == "concat_seqs" or in_args.command == "delete_metadata" or in_args.command == "pull_records_with_feature":
            timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.gb --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_sb[in_args.command]), number=int(in_args.iteration))
            sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
        elif in_args.command == "degenerate_sequence" or in_args.command == "find_CpG" or in_args.command == "find_orfs" or in_args.command == "find_restriction_sites":
            timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_nuc.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_sb[in_args.command]), number=int(in_args.iteration))
            sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
        elif in_args.command == "delete_features" or in_args.command == "extract_feature_sequences" or in_args.command == "extract_regions":
            timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.gb --%s \'%s\'", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_sb[in_args.command]), number=int(in_args.iteration))
            sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
        elif in_args.command == "delete_records":
            timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.fa --%s \'%s\'", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_sb[in_args.command]), number=int(in_args.iteration))
            sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
        elif in_args.command == "guess_alphabet" or in_args.command == "guess_format" or in_args.command == "map_features_nucl2prot" or in_args.command == "map_features_prot2nucl" or in_args.command == "merge":
            timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy \'%s\' --%s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (opts_sb[in_args.command], in_args.command), number=int(in_args.iteration))
            sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
        else:
            timer = timeit.timeit('from subprocess import Popen, PIPE; Popen("seqbuddy All_pannexins_pep.fa --%s %s", shell=True, stderr=PIPE, stdout=PIPE).communicate()' % (in_args.command, opts_sb[in_args.command]), number=int(in_args.iteration))
            sys.stdout.write("%s\n" % round(timer / int(in_args.iteration), 3))
    else:
        sys.stdout.write("Invalid command\n")
