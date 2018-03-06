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
from time import time
from subprocess import Popen, PIPE


class TempDir(object):
    def __init__(self):
        self.dir = next(self._make_dir())
        self.path = self.dir.name

    def _make_dir(self):
        temp_dir = TemporaryDirectory()
        yield temp_dir
        shutil.rmtree(self.path)


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


def pto(cmd, timeout, pipe=False):
    if pipe:
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    else:
        p = Popen(cmd, shell=True)
    start = round(time())
    end = start + timeout - 1
    while round(time()) <= end:
        if p.poll() is not None:
            return p.communicate()
    p.kill()
    return False


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog="performanceScanner", description="Check function time",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("reference", help="Specify input DNA sequences in genbank format")

    parser.add_argument("-t", "--tools", nargs="+", default=["all"],
                        help="Specify the module(s) or tool(s) to run")
    parser.add_argument("-s", "--skip", nargs="+", help="Specify module(s) or tool(s) that will be ignored")
    parser.add_argument("-3p", "--third_party", action='store_true', help="Include tools that use third party software")

    parser.add_argument("-i", "--iterations", action='store', default=10, help="Specify number of timeit replicates")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print out the result of each tool")
    parser.add_argument("-p", "--pause", action="store_true",
                        help="Stop execution until 'return' key pressed (only works in combination with -v)")
    parser.add_argument("-to", "--timeout", action='store', default=31536000, type=int, help="Set max execution time")
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

    res_dir = "{0}{1}results{1}{2}{1}".format(os.path.dirname(os.path.realpath(__file__)), os.path.sep, ref_name)
    if not os.path.isdir(res_dir):
        os.makedirs(res_dir)

    if not os.path.isfile("%s%s.gb" % (ref_dir, ref_name)):
        print(" -> Copying DNA file")
        shutil.copy(in_args.reference, "%s%s.gb" % (ref_dir, ref_name))

    if not os.path.isfile("%s%s_pep.gb" % (ref_dir, ref_name)):
        print(" -> Creating protein file")
        sb_pep = Sb.translate_cds(Sb.make_copy(seqbuddy), quiet=True)
        sb_pep.write("%s%s_pep.gb" % (ref_dir, ref_name))
        del sb_pep

    if not os.path.isfile("%s%s_rna.gb" % (ref_dir, ref_name)):
        print(" -> Creating RNA file")
        sb_rna = Sb.dna2rna(Sb.make_copy(seqbuddy))
        sb_rna.write("%s%s_rna.gb" % (ref_dir, ref_name))
        del sb_rna

    if not os.path.isfile("%s%s_aln.gb" % (ref_dir, ref_name)):
        print(" -> Creating alignment file")
        alignbuddy = Alb.faux_alignment(Sb.make_copy(seqbuddy), r_seed=12345)
        alignbuddy.write("%s%s_aln.gb" % (ref_dir, ref_name))
    else:
        alignbuddy = Alb.AlignBuddy("%s%s_aln.gb" % (ref_dir, ref_name))

    if not os.path.isfile("%s%s_pep_aln.gb" % (ref_dir, ref_name)):
        print(" -> Creating protein alignment file")
        alb_pep = Alb.faux_alignment(Sb.SeqBuddy("%s%s_pep.gb" % (ref_dir, ref_name)))
        alb_pep.write("%s%s_pep_aln.gb" % (ref_dir, ref_name))
        del alb_pep

    if not os.path.isfile("%s%s_rna_aln.gb" % (ref_dir, ref_name)):
        print(" -> Creating RNA alignment file")
        alb_rna = Alb.dna2rna(Alb.make_copy(alignbuddy))
        alb_rna.write("%s%s_rna_aln.gb" % (ref_dir, ref_name))
        del alb_rna

    if not os.path.isfile("%s%s_tree.nwk" % (ref_dir, ref_name)):
        print(" -> Creating tree file")
        from dendropy.simulate import treesim
        tree = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, num_extant_tips=len(seqbuddy))
        tree = tree.as_string("newick")
        for indx, rec in enumerate(seqbuddy.records):
            tree = re.sub("T%s:" % indx, "%s:" % rec.id, tree)
        phylobuddy = Pb.PhyloBuddy(tree)
        phylobuddy.write("%s%s_tree.nwk" % (ref_dir, ref_name))
        del tree
        del phylobuddy

    del seqbuddy
    del alignbuddy

    tmp_dir = TempDir()

    # Create all of the Tool objects for processing
    pd_tools = pd.read_csv("tools.csv", comment="#", escapechar="\\")
    tools = [Tool(tl.flag, tl.options, tl.module, ref_name, tl.reference, tl.third_party)
             for indx, tl in pd_tools.iterrows()]

    skip = [] if not in_args.skip else in_args.skip[0]

    # Benchmark each tool
    for tool in tools:
        if not in_args.third_party and tool.third_party:
            continue
        if tool.module in skip or tool.flag in skip:
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
                print("\033[92m%s (%s)\033[39m" % (tool, tool.module))
                verbose = ""
            else:
                sys.stdout.write("%s (%s): " % (tool.flag, tool.module))
                sys.stdout.flush()
                verbose = ", pipe=True"

            command = 'from performance_scanner import pto; '
            command += 'pto("{0} {1} --{2} {3} &> {4}{0}-{2}", '.format(tool.module, os.path.abspath(tool.reference),
                                                                        tool.flag, tool.options, res_dir)
            command += 'timeout=%s%s)' % (in_args.timeout, verbose)

            timer = timeit.timeit(command, number=int(in_args.iterations))
            sys.stdout.write("%s\n" % round(timer / int(in_args.iterations), 3))

            if in_args.verbose and in_args.pause:
                input("\033[91mPress 'return' to continue\033[39m")
