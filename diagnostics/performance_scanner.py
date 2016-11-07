#!/usr/bin/env python3
import sys
import os
import shutil
import argparse
import timeit
import re
import buddysuite.SeqBuddy as Sb
import buddysuite.AlignBuddy as Alb
import buddysuite.PhyloBuddy as Pb
from Bio.Alphabet import IUPAC
from tempfile import TemporaryDirectory
import pandas as pd


class TempDir(object):
    def __init__(self):
        self.dir = next(self._make_dir())
        self.path = self.dir.name

    def _make_dir(self):
        tmp_dir = TemporaryDirectory()
        yield tmp_dir
        rmtree(self.path)


class Tool(object):
    def __init__(self, flag, options, module, prefix, ref, third_party=False):
        self.flag = flag
        self.options = options if str(options) != "nan" else ""
        self.module = module
        self.reference = "reference%s%s" % (os.sep, prefix)
        if self.module == "phylobuddy":
            self.reference += "_tree.nwk"
        else:
            if ref == "pep":
                self.reference += "_pep"
            elif ref == "rna":
                self.reference += "_rna"
            elif ref == "dna/pep":
                self.reference += ".gb %s_pep" % self.reference
            elif ref == "dna/dna":
                self.reference += ".gb %s" % self.reference
            if self.module == "alignbuddy":
                self.reference += "_aln"
            self.reference += ".gb"
        self.third_party = third_party

    def __str__(self):
        return "$: %s %s --%s %s" % (self.module, self.reference, self.flag, self.options)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="performanceScanner", description="Check function time",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("reference", help="Specify input DNA sequences in genbank format")

    parser.add_argument("-t", "--tools", nargs="+", default=["all"],
                        help="Specify the module(s) or tool(s) to run")
    parser.add_argument("-3p", "--third_party", action='store_true', help="Include tools that use third party software")

    parser.add_argument("-i", "--iterations", action='store', default=10, help="Specify number of timeit replicates")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print out the result of each tool")
    parser.add_argument("-p", "--pause", action="store_true",
                        help="Stop execution until 'return' key pressed (only workes in combination with -v)")
    in_args = parser.parse_args()

    # Validate input reference file
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

    # Create or load all necessary reference files
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

    tmp_dir = TempDir()

    # Create all of the Tool objects for processing
    pd_tools = pd.read_csv("tools.csv", comment="#", escapechar="\\")
    tools = [Tool(tl.flag, tl.options, tl.module, ref_name, tl.reference, tl.third_party) for indx, tl in pd_tools.iterrows()]

    # Benchmark each tool
    for tool in tools:
        if not in_args.third_party and tool.third_party:
            continue

        if any(i in in_args.tools for i in ["all", tool.flag, tool.module]):  # Allows multiple tools to be called
            # Catch any hooks in the options and convert
            hook = re.search("__(.+)__", tool.options)
            if hook:
                hook = hook.group(1)
                if "ref" in hook:
                    tool.options = re.sub("__.+__", "reference%s%s%s" % (os.sep, ref_name, hook[3:]), tool.options)
                elif hook == "tmp":
                    tool.options = re.sub("__.+__", tmp_dir.path, tool.options)

            if in_args.verbose:
                print(tool)
                pipe = ""
            else:
                sys.stdout.write("%s: " % tool.flag)
                sys.stdout.flush()
                pipe = ", stderr=PIPE, stdout=PIPE"

            command = 'from subprocess import Popen, PIPE; '
            command += 'Popen("%s %s --%s %s", ' % (tool.module, tool.reference, tool.flag, tool.options)
            command += 'shell=True%s).communicate()' % pipe

            timer = timeit.timeit(command, number=int(in_args.iterations))
            sys.stdout.write("%s\n" % round(timer / int(in_args.iterations), 3))

            if in_args.verbose and in_args.pause:
                input("Press 'return' to continue")

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