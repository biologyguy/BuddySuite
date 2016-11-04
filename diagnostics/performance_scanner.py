#!/usr/bin/env python3
import sys
import os
import shutil
import argparse
import timeit
import math
import buddysuite.SeqBuddy as Sb
import buddysuite.AlignBuddy as Alb
import buddysuite.PhyloBuddy as Pb
from Bio.Alphabet import IUPAC
from tempfile import TemporaryDirectory


class TempDir(object):
    def __init__(self):
        self.dir = next(self._make_dir())
        self.path = self.dir.name

    def _make_dir(self):
        tmp_dir = TemporaryDirectory()
        yield tmp_dir
        rmtree(self.path)


class Tool(object):
    def __init__(self, flag, options, module, ref, third_party=False):
        self.flag = flag
        self.options = options
        self.module = module
        self.reference = ref
        self.third_party = third_party

    def ref_file(self, file_prefix):
        output = file_prefix

        if self.module == "phylobuddy":
            return output + "_tree.nwk"

        if self.reference == "pep":
            output += "_pep"
        elif self.reference == "rna":
            output += "_rna"

        if self.module == "alignbuddy":
            output += "_aln"

        return output + ".gb"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="performanceScanner", description="Check function time", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("reference", help="Specify input DNA sequences in genbank format")

    parser.add_argument("-t", "--tools", help="Specify the module(s) or tool(s) to run", action='store', default="all")
    parser.add_argument("-3p", "--third_party", help="Include tools that use third party software", action='store_true')

    parser.add_argument("-i", "--iterations", help="Specify number of timeit replicates", action='store', default=10)
    in_args = parser.parse_args()

    if not os.path.isfile(in_args.reference):
        sys.stderr("Error: Reference file does not exist\n")
        sys.exit()

    seqbuddy = Sb.SeqBuddy(in_args.reference)
    if seqbuddy.alpha != IUPAC.ambiguous_dna:
        sys.stderr("Error: Reference file must be DNA\n")
        sys.exit()

    if seqbuddy.in_format not in ["genbank", "gb"]:
        sys.stderr("Error: Reference file must be GenBank format\n")
        sys.exit()

    tmp_dir = TempDir()
    ref_dir = "{0}{1}reference{1}".format(os.path.dirname(os.path.realpath(__file__)), os.path.sep)
    ref_name = in_args.reference.split(os.sep)[-1]
    ref_name = os.path.splitext(ref_name)[0]

    if not os.path.isfile("%s%s.gb" % (ref_dir, ref_name)):
        shutil.copy(in_args.reference, "%s%s.gb" % (ref_dir, ref_name))

    if not os.path.isfile("%s%s_pep.gb" % (ref_dir, ref_name)):
        sb_pep = Sb.translate_cds(Sb.make_copy(seqbuddy))
        sb_pep.write("%s%s_pep.gb" % (ref_dir, ref_name))

    if not os.path.isfile("%s%s_rna.gb" % (ref_dir, ref_name)):
        sb_rna = Sb.dna2rna(Sb.make_copy(seqbuddy))
        sb_rna.write("%s%s_rna.gb" % (ref_dir, ref_name))

    if not os.path.isfile("%s%s_aln.gb" % (ref_dir, ref_name)):
        alignbuddy = Alb.generate_msa(Sb.make_copy(seqbuddy), "mafft")
        alignbuddy.write("%s%s_aln.gb" % (ref_dir, ref_name))
    else:
        alignbuddy = Alb.AlignBuddy("%s%s_aln.gb" % (ref_dir, ref_name))

    if not os.path.isfile("%s%s_pep_aln.gb" % (ref_dir, ref_name)):
        alb_pep = Alb.translate_cds(Alb.make_copy(alignbuddy))
        alb_pep.write("%s%s_pep_aln.gb" % (ref_dir, ref_name))

    if not os.path.isfile("%s%s_rna_aln.gb" % (ref_dir, ref_name)):
        alb_rna = Alb.dna2rna(Alb.make_copy(alignbuddy))
        alb_rna.write("%s%s_rna_aln.gb" % (ref_dir, ref_name))

    if not os.path.isfile("%s%s_tree.nwk" % (ref_dir, ref_name)):
        phylobuddy = Pb.generate_tree(Alb.make_copy(alignbuddy), "fasttree")
        phylobuddy.write("%s%s_tree.nwk" % (ref_dir, ref_name))
    else:
        phylobuddy = Pb.PhyloBuddy("%s%s_tree.nwk" % (ref_dir, ref_name))

    tools = []
    tools.append(Tool("annotate", "misc_feature 1-10", "seqbuddy", "dna", False))
    tools.append(Tool("ave_seq_length", "", "seqbuddy", "dna", False))
    tools.append(Tool("back_translate", "", "seqbuddy", "pep", False))
    tools.append(Tool("bl2seq", "", "seqbuddy", "dna", True))
    tools.append(Tool("blast", "mnemiopsis_leidyi", "seqbuddy", "dna", True))
    tools.append(Tool("clean_seq", "", "seqbuddy", "dna", False))
    tools.append(Tool("complement", "", "seqbuddy", "dna", False))
    tools.append(Tool("concat_seqs", "", "seqbuddy", "dna", False))
    tools.append(Tool("count_codons", "", "seqbuddy", "dna", False))
    tools.append(Tool("count_residues", "concatenate", "seqbuddy", "dna", False))
    tools.append(Tool("degenerate_sequence", "1", "seqbuddy", "pep", False))
    tools.append(Tool("delete_features", "TMD1", "seqbuddy", "dna", False))
    tools.append(Tool("delete_large", "4000", "seqbuddy", "dna", False))
    tools.append(Tool("delete_metadata", "", "seqbuddy", "dna", False))
    tools.append(Tool("delete_records", "Aae-", "seqbuddy", "dna", False))
    tools.append(Tool("delete_repeats", "", "seqbuddy", "dna", False))
    tools.append(Tool("delete_small", "1000", "seqbuddy", "dna", False))
    tools.append(Tool("extract_feature_sequences", "TMD1", "seqbuddy", "dna", False))
    tools.append(Tool("extract_regions", ",50", "seqbuddy", "dna", False))
    tools.append(Tool("find_CpG", "-o genbank", "seqbuddy", "dna", False))
    tools.append(Tool("find_orfs", "-o genbank", "seqbuddy", "dna", False))
    tools.append(Tool("find_pattern", "\'LLL\' \'LLY\'", "seqbuddy", "pep", False))
    tools.append(Tool("find_repeats", "3", "seqbuddy", "dna", False))
    tools.append(Tool("find_restriction_sites", "commercial", "seqbuddy", "dna", False))
    tools.append(Tool("group_by_prefix",  tmp_dir.path, "seqbuddy", "dna", False))
    tools.append(Tool("group_by_regex", "%s Sra" % tmp_dir.path, "seqbuddy", "dna", False))
    tools.append(Tool("guess_alphabet", "reference/*", "seqbuddy", "dna", False))
    tools.append(Tool("guess_format", "reference/*", "seqbuddy", "dna", False))
    tools.append(Tool("hash_seq_ids", "10", "seqbuddy", "dna", False))
    tools.append(Tool("insert_seq", "DYKDDDDK 0", "seqbuddy", "pep", False))
    tools.append(Tool("isoelectric_point", "", "seqbuddy", "pep", False))
    tools.append(Tool("list_features", "", "seqbuddy", "dna", False))
    tools.append(Tool("list_ids", "3", "seqbuddy", "dna", False))
    tools.append(Tool("lowercase", "", "seqbuddy", "dna", False))
    tools.append(Tool("make_ids_unique", "10 \'-\'", "seqbuddy", "dna", False))
    tools.append(Tool("map_features_nucl2prot", "", "seqbuddy", "dna/pep", False))
    tools.append(Tool("map_features_prot2nucl", "", "seqbuddy", "dna/pep", False))
    tools.append(Tool("merge", "", "seqbuddy", "dna/dna", False))
    tools.append(Tool("molecular_weight", "", "seqbuddy", "dna", False))
    tools.append(Tool("num_seqs", "", "seqbuddy", "dna", False))
    tools.append(Tool("order_features_alphabetically", "", "seqbuddy", "dna", False))
    tools.append(Tool("order_features_by_position", "", "seqbuddy", "dna", False))
    tools.append(Tool("order_ids", "", "seqbuddy", "dna", False))
    tools.append(Tool("order_ids_randomly", "", "seqbuddy", "dna", False))
    tools.append(Tool("prosite_scan", "strict", "seqbuddy", "pep", True))
    tools.append(Tool("pull_random_record", math.ceil(len(seqbuddy) / 10), "seqbuddy", "dna", False))
    tools.append(Tool("pull_records", "Aae-Panx", "seqbuddy", "dna", False))
    tools.append(Tool("pull_record_ends", "-100", "seqbuddy", "dna", False))
    tools.append(Tool("pull_records_with_feature", "PHOSPHO", "seqbuddy", "dna", False))
    tools.append(Tool("purge", "230", "seqbuddy", "dna", True))
    tools.append(Tool("rename_ids", "Mle Mnemiopsis", "seqbuddy", "dna", False))
    tools.append(Tool("replace_subseq", "IL].{1,4}[IL] _motif_", "seqbuddy", "pep", False))
    tools.append(Tool("reverse_complement", "", "seqbuddy", "dna", False))
    tools.append(Tool("reverse_transcribe", "", "seqbuddy", "rna", False))
    tools.append(Tool("screw_formats", "fasta", "seqbuddy", "dna", False))
    tools.append(Tool("select_frame", "3", "seqbuddy", "dna", False))
    tools.append(Tool("shuffle_seqs", "", "seqbuddy", "dna", False))
    tools.append(Tool("translate", "", "seqbuddy", "dna", False))
    tools.append(Tool("translate6frames", "", "seqbuddy", "dna", False))
    tools.append(Tool("transcribe", "", "seqbuddy", "dna", False))
    tools.append(Tool("transmembrane_domains", "", "seqbuddy", "pep", True))
    tools.append(Tool("uppercase", "", "seqbuddy", "dna", False))
    
    tools.append(Tool("alignment_lengths", "", "alignbuddy", "", False))
    tools.append(Tool("back_transcribe", "", "alignbuddy", "", False))
    tools.append(Tool("bootstrap", "3", "alignbuddy", "", False))
    tools.append(Tool("clean_seq", "strict", "alignbuddy", "", False))
    tools.append(Tool("concat_alignments", "[a-z]{3}-Panx", "alignbuddy", "", False))
    tools.append(Tool("consensus", "", "alignbuddy", "", False))
    tools.append(Tool("delete_records", "PanxγA", "alignbuddy", "", False))
    tools.append(Tool("enforce_triplets", "", "alignbuddy", "", False))
    tools.append(Tool("extract_regions", "10 110", "alignbuddy", "", False))
    tools.append(Tool("generate_alignment", "mafft", "alignbuddy", "", False))
    tools.append(Tool("hash_ids", "10", "alignbuddy", "", False))
    tools.append(Tool("list_ids", "3", "alignbuddy", "", False))
    tools.append(Tool("lowercase", "", "alignbuddy", "", False))
    tools.append(Tool("mapfeat2align", "All_pannexins_pep.gb", "alignbuddy", "", False))
    tools.append(Tool("num_seqs", "", "alignbuddy", "", False))
    tools.append(Tool("order_ids", "rev", "alignbuddy", "", False))
    tools.append(Tool("pull_records", "PanxβB PanxβC", "alignbuddy", "", False))
    tools.append(Tool("rename_ids", "Mle Mnemiopsis", "alignbuddy", "", False))
    tools.append(Tool("screw_formats", "fasta", "alignbuddy", "", False))
    tools.append(Tool("split_to_files", "~/BuddySuite/diagnostics/files_al Pannexin_", "alignbuddy", "", False))
    tools.append(Tool("transcribe", "", "alignbuddy", "", False))
    tools.append(Tool("translate", "", "alignbuddy", "", False))
    tools.append(Tool("trimal", "clean", "alignbuddy", "", False))
    tools.append(Tool("uppercase", "", "alignbuddy", "", False))
    
    tools.append(Tool("collapse_polytomies", "length 0.1", "phylobuddy", "", False))
    tools.append(Tool("consensus_tree", "0.5", "phylobuddy", "", False))
    tools.append(Tool("display_trees", "", "phylobuddy", "", False))
    tools.append(Tool("distance", "uwrf", "phylobuddy", "", False))
    tools.append(Tool("generate_tree", "raxmlHPC-SSE3 -o nexus", "phylobuddy", "", False))
    tools.append(Tool("hash_ids", "3", "phylobuddy", "", False))
    tools.append(Tool("list_ids", "3", "phylobuddy", "", False))
    tools.append(Tool("num_tips", "", "phylobuddy", "", False))
    tools.append(Tool("print_trees", "", "phylobuddy", "", False))
    tools.append(Tool("prune_taxa", "Cte-PanxζI", "phylobuddy", "", False))
    tools.append(Tool("rename_ids", "α A", "phylobuddy", "", False))
    tools.append(Tool("root", "2", "phylobuddy", "", False))
    tools.append(Tool("screw_formats", "newick", "phylobuddy", "", False))
    tools.append(Tool("show_unique", "", "phylobuddy", "", False))
    tools.append(Tool("split_polytomies", "", "phylobuddy", "", False))
    tools.append(Tool("unroot", "", "phylobuddy", "", False))

    for tool in tools[:5]:
        assert tool.reference in ["dna", "dna/dna", "pep", "dna/pep", "rna", "tree"]

        if in_args.tools in ["all", tool.flag, tool.module]:
            command = 'from subprocess import Popen, PIPE; '
            command += 'Popen("%s %s --%s %s", ' % (tool.module, tool.ref_file(ref_name), tool.flag, tool.options)
            command += 'shell=True, stderr=PIPE, stdout=PIPE).communicate()'

            sys.stdout.write("%s: " % tool.flag)
            sys.stdout.flush()
            timer = timeit.timeit(command, number=int(in_args.iterations))
            sys.stdout.write("%s\n" % round(timer / int(in_args.iterations), 3))

    sys.exit()
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