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
        temp_dir = TemporaryDirectory()
        yield temp_dir
        rmtree(self.path)


class Tool(object):
    def __init__(self, flag, options, module, prefix, ref, third_party=False):
        self.flag = flag
        self.options = options if str(options) != "nan" else ""
        self.module = module
        self.reference = "reference%s%s" % (os.sep, prefix)
        if self.module == "phylobuddy" and "tree" in ref:
            self.reference += "_tree.nwk"
            if ref == "tree/tree":
                self.reference += " %s" % self.reference
        else:
            if ref == "pep":
                self.reference += "_pep"
            elif ref == "rna":
                self.reference += "_rna"
            elif ref == "dna/pep":
                self.reference += ".gb %s_pep" % self.reference
            elif ref == "dna/dna":
                self.reference += ".gb %s" % self.reference
            self.reference += ".gb"

            if self.module in ["alignbuddy", "phylobuddy"] and flag != "generate_alignment":
                self.reference = re.sub("\.gb", "_aln.gb", self.reference)

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
                print("\033[92m%s\033[39m" % tool)
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
                input("\033[91mPress 'return' to continue\033[39m")
