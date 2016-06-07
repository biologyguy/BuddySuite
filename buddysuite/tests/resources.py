# coding=utf-8
import sys
import os
from copy import deepcopy

main_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, main_dir)
sys.path.insert(0, os.path.abspath(".."))

try:
    from buddysuite import AlignBuddy as Alb
    from buddysuite import SeqBuddy as Sb
    from buddysuite import PhyloBuddy as Pb
    from buddysuite import DatabaseBuddy as Db
    from buddysuite import buddy_resources as br
    from buddysuite import MyFuncs
except ImportError:
    import AlignBuddy as Alb
    import SeqBuddy as Sb
    import PhyloBuddy as Pb
    import DatabaseBuddy as Db
    import buddy_resources as br
    import MyFuncs

# This file (conftest.py) must be in the same directory as unit_test_resources
RESOURCE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'unit_test_resources')


class Resources(object):
    """
    Resources are organized by molecule, number of alignmentts, and file format

    self.resource_list[<molecule_type>][<file_format>][<num_aligns>]
    <molecule_type>:
        'dna', 'rna', or 'pep'

    <num_aligns>:
        'multi' or 'single'

    <file_format>:
        'clustal'
        'fasta'
        'gb'
        'nexus'
        'phylip'
        'phylipr'
        'phylipss'
        'phylipsr'
        'stockholm'
    """

    def __init__(self):
        base_dict_structure = {'dna': {'single': {}, 'multi': {}},
                               'rna': {'single': {}, 'multi': {}},
                               'pep': {'single': {}, 'multi': {}}}

        self.resources = deepcopy(base_dict_structure)
        self.resources['dna']['single'] = {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
            ("clustal", "{path}/Mnemiopsis_cds.clus"),
            ("fasta", "{path}/Mnemiopsis_cds_aln.fa"),
            ("gb", "{path}/Mnemiopsis_cds_aln.gb"),
            ("nexus", "{path}/Mnemiopsis_cds.nex"),
            ("phylip", "{path}/Mnemiopsis_cds.phy"),
            ("phylipr", "{path}/Mnemiopsis_cds.phyr"),
            ("phylipss", "{path}/Mnemiopsis_cds.physs"),
            ("phylipsr", "{path}/Mnemiopsis_cds.physr"),
            ("stockholm", "{path}/Mnemiopsis_cds.stklm")]}

        self.resources['dna']['multi'] = {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
            ("clustal", "{path}/Alignments_cds.clus"),
            ("phylip", "{path}/Alignments_cds.phy"),
            ("phylipr", "{path}/Alignments_cds.phyr"),
            ("phylipss", "{path}/Alignments_cds.physs"),
            ("phylipsr", "{path}/Alignments_cds.physr"),
            ("stockholm", "{path}/Alignments_cds.stklm")]}
        self.resources['rna']['single'] = {"nexus": "{path}/Mnemiopsis_rna.nex".format(path=RESOURCE_PATH)}
        self.resources['pep']['single'] = {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
            ("gb", "{path}/Mnemiopsis_pep_aln.gb"),
            ("nexus", "{path}/Mnemiopsis_pep.nex"),
            ("phylip", "{path}/Mnemiopsis_pep.phy"),
            ("phylipr", "{path}/Mnemiopsis_pep.phyr"),
            ("phylipss", "{path}/Mnemiopsis_pep.physs"),
            ("phylipsr", "{path}/Mnemiopsis_pep.physr"),
            ("stockholm", "{path}/Mnemiopsis_pep.stklm")]}
        self.resources['pep']['multi'] = {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
            ("clustal", "{path}/Alignments_pep.clus"),
            ("phylip", "{path}/Alignments_pep.phy"),
            ("phylipr", "{path}/Alignments_pep.phyr"),
            ("phylipss", "{path}/Alignments_pep.physs"),
            ("phylipsr", "{path}/Alignments_pep.physr"),
            ("stockholm", "{path}/Alignments_pep.stklm")]}

        # Create new AlignBuddy objects for each resource file
        self.alb_objs = deepcopy(base_dict_structure)
        for mol in self.resources:
            for num in self.resources[mol]:
                for file_format in self.resources[mol][num]:
                    self.alb_objs[mol][num][file_format] = Alb.AlignBuddy(self.resources[mol][num][file_format])

        self.code_dict = {"molecule": {"p": "pep", "d": "dna", "r": "rna"},
                          "num_aligns": {"o": "single", "m": "multi"},
                          "format": {"c": "clustal", "f": "fasta", "g": "gb", "n": "nexus", "py": "phylip",
                                     "pr": "phylipr", "pss": "phylipss", "psr": "phylipsr", "s": "stockholm"}}

        self.single_letter_codes = {"p": "pep", "d": "dna", "r": "rna",
                                    "o": "single", "m": "multi",
                                    "c": "clustal", "f": "fasta", "g": "gb", "n": "nexus", "py": "phylip",
                                    "pr": "phylipr", "pss": "phylipss", "psr": "phylipsr", "s": "stockholm"}

    def parse_code(self, code="", strict=False):
        """
        Take in the letter codes for a query and determine the final groups to be returned
        When codes from a particular category are ommited, pull in all possibilities for that categroy
        :param code: Letter codes (explained in Class definition)
        :type code: str
        :param strict: Only return the exact resources included in code (don't inflate empty types)
        :type strict: bool
        :return: The complete group of resources to be used
        :rtype: dict
        """
        results = {"molecule": [], "num_aligns": [], "format": []}
        code = code.split()
        # Sorry about this maddness.. Each code is checked against each of the types in self.code_dict
        # and pushed into the final results if it is found there.
        for i in code:
            for j in results:
                if i in self.code_dict[j]:
                    results[j].append(i)

        # Fill up fields with all possibilities if nothing is given
        for result_type in results:
            if not results[result_type] and not strict:
                results[result_type] = [key for key in self.code_dict[result_type]]
        return results

    def get_key(self, code):
        code = code.split()
        if len(code) != 3:
            raise AttributeError("Only explicit three-component codes are accepted")
        for letter in code:
            if letter not in self.single_letter_codes:
                raise AttributeError("Malformed letter code, '%s' not recognized" % letter)

        output = {}
        for component in ["molecule", "num_aligns", "format"]:
            for indx, letter in enumerate(code):
                if letter in self.code_dict[component]:
                    if component in output:
                        raise AttributeError("Malformed letter code, trying to append multiple values to the %s "
                                             "component of the key" % component)
                    output[component] = self.single_letter_codes[letter]
                    del code[indx]
                    break

        # return [output["molecule"], output["num_aligns"], output["format"]]
        return output

    def get(self, code="", mode="objs"):
        """
        Returns copies of AlignBuddy objects of the path to their resource files
        :param code: Letter codes (explained in Class definition)
        :type code: str
        :param mode: Return either AlignBuddy "objs" (default) or "paths"
        :type mode: str
        :return: AlignBuddy objects or resource paths as controlled by mode {key: resource}
        :rtype: dict
        """
        files = self.parse_code(code)
        output = {}
        slc = self.single_letter_codes
        for molecule in files["molecule"]:
            for num_aligns in files["num_aligns"]:
                for _format in files["format"]:
                    try:
                        if mode == "paths":
                            new_obj = self.resources[slc[molecule]][slc[num_aligns]][slc[_format]]
                        elif mode == "objs":
                            new_obj = self.alb_objs[slc[molecule]][slc[num_aligns]][slc[_format]]
                            new_obj = Alb.make_copy(new_obj)
                        else:
                            raise ValueError("The 'mode' parameter only accepts 'objs' or 'paths' as input.")
                        output["%s %s %s" % (num_aligns, molecule, _format)] = new_obj
                    except KeyError:
                        pass
        return output

    def get_list(self, code="", mode="objs"):
        return [value for key, value in self.get(code=code, mode=mode).items()]

    def get_one(self, code, mode="objs"):
        if len(code.split()) != 3:
            raise AttributeError("Only explicit three-component codes are accepted")
        output = self.get_list(code, mode)
        return None if not output else output[0]
