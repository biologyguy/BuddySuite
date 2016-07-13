# coding=utf-8
import sys
import os
from copy import deepcopy
import re
from hashlib import md5

from .. import AlignBuddy as Alb
from .. import SeqBuddy as Sb
from .. import PhyloBuddy as Pb
from .. import DatabaseBuddy as Db
from .. import buddy_resources as br

# This file (conftest.py) must be in the same directory as unit_test_resources
RESOURCE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'unit_test_resources')


# #################################  -  SeqBuddy  -  ################################## #
class SbResources(object):
    def __init__(self):
        base_dict_structure = {'dna': {}, 'rna': {}, 'pep': {}}

        self.resources = deepcopy(base_dict_structure)
        self.resources['dna'] = {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
            ("clustal", "{path}/Mnemiopsis_cds.clus"),
            ("fasta", "{path}/Mnemiopsis_cds.fa"),
            ("gb", "{path}/Mnemiopsis_cds.gb"),
            ("nexus", "{path}/Mnemiopsis_cds.nex"),
            ("phylip", "{path}/Mnemiopsis_cds.phy"),
            ("phylipr", "{path}/Mnemiopsis_cds.phyr"),
            ("phylipss", "{path}/Mnemiopsis_cds.physs"),
            ("phylipsr", "{path}/Mnemiopsis_cds.physr"),
            ("stockholm", "{path}/Mnemiopsis_cds.stklm"),
            ("embl", "{path}/Mnemiopsis_cds.embl"),
            ("seqxml", "{path}/Mnemiopsis_cds.seqxml")]}

        self.resources['rna'] = {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
            ("fasta", "{path}/Mnemiopsis_rna.fa"),
            ("nexus", "{path}/Mnemiopsis_rna.nex")]}
        self.resources['pep'] = {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
            ("fasta", "{path}/Mnemiopsis_pep.fa"),
            ("gb", "{path}/Mnemiopsis_pep.gb"),
            ("nexus", "{path}/Mnemiopsis_pep.nex"),
            ("phylip", "{path}/Mnemiopsis_pep.phy"),
            ("phylipr", "{path}/Mnemiopsis_pep.phyr"),
            ("phylipss", "{path}/Mnemiopsis_pep.physs"),
            ("phylipsr", "{path}/Mnemiopsis_pep.physr"),
            ("stockholm", "{path}/Mnemiopsis_pep.stklm")]}

        # Create new SeqBuddy objects for each resource file
        self.sb_objs = deepcopy(base_dict_structure)
        for mol in self.resources:
            for file_format in self.resources[mol]:
                self.sb_objs[mol][file_format] = Sb.SeqBuddy(self.resources[mol][file_format])

        self.code_dict = {"molecule": {"p": "pep", "d": "dna", "r": "rna"},
                          "format": {"c": "clustal", "f": "fasta", "g": "gb", "n": "nexus", "py": "phylip",
                                     "pr": "phylipr", "pss": "phylipss", "psr": "phylipsr", "s": "stockholm",
                                     "e": "embl", "x": "seqxml"}}

        self.single_letter_codes = {"p": "pep", "d": "dna", "r": "rna",
                                    "c": "clustal", "f": "fasta", "g": "gb", "n": "nexus", "py": "phylip",
                                    "pr": "phylipr", "pss": "phylipss", "psr": "phylipsr", "s": "stockholm",
                                    "e": "embl", "x": "seqxml"}
        self.res_path = RESOURCE_PATH

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
        results = {"molecule": [], "format": []}
        code = code.split()
        for i in code:
            for j in results:
                if i in self.code_dict[j]:
                    results[j].append(i)

        # Fill up a field with all possibilities if nothing is given
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
        for component in ["molecule", "format"]:
            for indx, letter in enumerate(code):
                if letter in self.code_dict[component]:
                    if component in output:
                        raise AttributeError("Malformed letter code, trying to append multiple values to the %s "
                                             "component of the key" % component)
                    output[component] = self.single_letter_codes[letter]
                    del code[indx]
                    break

        # return [output["molecule"], output["format"]]
        return output

    def get(self, code="", mode="objs"):
        """
        Returns copies of SeqBuddy objects
        :param code:
        :param mode: {"objs", "paths"}
        :return: OrderedDict {key: resource}
        """
        files = self.parse_code(code)
        output = {}
        slc = self.single_letter_codes
        for molecule in files["molecule"]:
            for _format in files["format"]:
                try:
                    if mode == "paths":
                        new_obj = self.resources[slc[molecule]][slc[_format]]
                    elif mode == "objs":
                        new_obj = self.sb_objs[slc[molecule]][slc[_format]]
                        new_obj = Sb.make_copy(new_obj)
                    else:
                        raise ValueError("The 'mode' parameter only accepts 'objs' or 'paths' as input.")
                    output["%s %s" % (molecule, _format)] = new_obj
                except KeyError:
                    pass
        return output

    def get_list(self, code="", mode="objs"):
        return [value for key, value in self.get(code=code, mode=mode).items()]

    def get_one(self, code, mode="objs"):
        output = self.get_list(code, mode)
        return None if not output or len(output) > 1 else output[0]

    def deets(self, code):
        code = code.split()
        return {"type": self.code_dict["type"][code[0]],
                "format": br.parse_format(self.code_dict["format"][code[1]])}


class SbHelpers(object):
    def __init__(self):
        self.resource_path = RESOURCE_PATH
        self.write_file = br.TempFile()

    def seqs2hash(self, _seqbuddy, mode='hash'):
        if _seqbuddy.out_format in ["gb", "genbank"]:
            for _rec in _seqbuddy.records:
                try:
                    if re.search("(\. )+", _rec.annotations['organism']):
                        _rec.annotations['organism'] = "."
                except KeyError:
                    pass

        if _seqbuddy.out_format == "phylipsr":
            self.write_file.write(br.phylip_sequential_out(_seqbuddy, relaxed=True, _type="seqbuddy"))
        elif _seqbuddy.out_format == "phylipss":
            self.write_file.write(br.phylip_sequential_out(_seqbuddy, relaxed=False, _type="seqbuddy"))
        else:
            _seqbuddy.write(self.write_file.path)

        seqs_string = "{0}\n".format(self.write_file.read().rstrip())
        self.write_file.clear()

        if mode != "hash":
            return seqs_string

        _hash = md5(seqs_string.encode()).hexdigest()
        return _hash

    @staticmethod
    def string2hash(_input):
        return md5(_input.encode("utf-8")).hexdigest()

    def features2hash(self, _seqbuddy, mode='hash'):
        feature_list = [f.id + ": " + str(f.features) for f in _seqbuddy.records]
        feature_string = "\n".join(feature_list)

        if mode != 'hash':
            return feature_string

        _hash = self.string2hash(feature_string)
        return _hash


# #################################  -  AlignBuddy  -  ################################ #
class AlbResources(object):
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
        self.res_path = RESOURCE_PATH

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


class AlbHelpers(object):
    def __init__(self):
        self.resource_path = RESOURCE_PATH

    @staticmethod
    def align2hash(alignbuddy=None, mode='hash'):
        if not alignbuddy:
            raise AttributeError("AlignBuddy object required")

        if mode != "hash":
            return "{0}".format(str(alignbuddy))
        _hash = md5("{0}".format(str(alignbuddy)).encode("utf-8")).hexdigest()
        return _hash

    @staticmethod
    def string2hash(_input):
        return md5(_input.encode("utf-8")).hexdigest()


# ################################  -  PhyloBuddy  -  ################################# #
class PbResources(object):
    def __init__(self):
        base_dict_structure = {'single': {}, 'multi': {}}

        self.resources = deepcopy(base_dict_structure)
        self.resources['single'] = {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
            ("newick", "{path}/single_tree.newick"),
            ("nexus", "{path}/single_tree.nex"),
            ("nexml", "{path}/single_tree.xml")]}

        self.resources['multi'] = {file_format: name.format(path=RESOURCE_PATH) for file_format, name in [
            ("newick", "{path}/multi_tree.newick"),
            ("nexus", "{path}/multi_tree.nex"),
            ("nexml", "{path}/multi_tree.xml")]}

        # Create new PhyloBuddy objects for each resrouce file
        self.pb_objs = deepcopy(base_dict_structure)
        for num in self.resources:
            for file_format in self.resources[num]:
                self.pb_objs[num][file_format] = Pb.PhyloBuddy(self.resources[num][file_format])

        self.code_dict = {"num_trees": {"o": "single", "m": "multi"},
                          "format": {"k": "newick", "n": "nexus", "l": "nexml"}}

        self.single_letter_codes = {"o": "single", "m": "multi",
                                    "k": "newick", "n": "nexus", "l": "nexml"}
        self.res_path = RESOURCE_PATH

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
        results = {"num_trees": [], "format": []}
        code = code.split()
        for i in code:
            for j in results:
                if i in self.code_dict[j]:
                    results[j].append(i)

        # Fill up a field with all possibilities if nothing is given
        for result_type in results:
            if not results[result_type] and not strict:
                results[result_type] = [key for key in self.code_dict[result_type]]
        return results

    def get_key(self, code):
        code = code.split()
        if len(code) != 2:
            raise AttributeError("Only explicit two-component codes are accepted")
        for letter in code:
            if letter not in self.single_letter_codes:
                raise AttributeError("Malformed letter code, '%s' not recognized" % letter)

        output = {}
        for component in ["num_trees", "format"]:
            for indx, letter in enumerate(code):
                if letter in self.code_dict[component]:
                    if component in output:
                        raise AttributeError("Malformed letter code, trying to append multiple values to the %s "
                                             "component of the key" % component)
                    output[component] = self.single_letter_codes[letter]
                    del code[indx]
                    break

        # return [output["num_trees"], output["format"]]
        return output

    def get(self, code="", mode="objs"):
        """
        Returns copies of PhyloBuddy objects, the
        :param code:
        :param mode: {"objs", "paths"}
        :return: OrderedDict {key: resource}
        """
        files = self.parse_code(code)
        output = {}
        slc = self.single_letter_codes
        for num_trees in files["num_trees"]:
            for _format in files["format"]:
                try:
                    if mode == "paths":
                        new_obj = self.resources[slc[num_trees]][slc[_format]]
                    elif mode == "objs":
                        new_obj = self.pb_objs[slc[num_trees]][slc[_format]]
                        new_obj = Pb.make_copy(new_obj)
                    else:
                        raise ValueError("The 'mode' parameter only accepts 'objs' or 'paths' as input.")
                    output["%s %s" % (num_trees, _format)] = new_obj
                except KeyError:
                    pass
        return output

    def get_list(self, code="", mode="objs"):
        return [value for key, value in self.get(code=code, mode=mode).items()]

    def get_one(self, code, mode="objs"):
        if len(code.split()) != 2:
            raise AttributeError("Only explicit two-component codes are accepted")
        output = self.get_list(code, mode)
        return None if not output else output[0]

    def deets(self, code):
        code = code.split()
        return {"num_trees": self.code_dict["num_trees"][code[0]],
                "format": br.parse_format(self.code_dict["format"][code[1]])}


class PbHelpers(object):
    def __init__(self):
        self.resource_path = RESOURCE_PATH
        self.write_file = br.TempFile()

    @staticmethod
    def phylo2hash(_phylobuddy, mode='hash'):
        if mode != "hash":
            return "{0}\n".format(str(_phylobuddy).rstrip())
        _hash = md5("{0}\n".format(str(_phylobuddy).rstrip()).encode('utf-8')).hexdigest()
        return _hash

    @staticmethod
    def string2hash(_input):
        return md5(_input.encode("utf-8")).hexdigest()

