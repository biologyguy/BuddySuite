#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from functools import partial
import shutil
from shutil import which, rmtree, copytree
from urllib import request, error
from platform import *
from os import path
from configparser import *
import copy
from tempfile import TemporaryDirectory
from subprocess import Popen
import zipfile
from inspect import getsourcefile
from collections import OrderedDict
import buddy_resources as br
import string
from random import choice

import argparse

_version = br.Version("BuddySuite", 1, 'alpha', br.contributors)

fmt = lambda prog: br.CustomHelpFormatter(prog)

parser = argparse.ArgumentParser(prog="buddysuite", formatter_class=fmt, add_help=False, usage=argparse.SUPPRESS,
                                 description='''\
\033[1mBuddySuite Installer\033[m
  Install, upgrade, or uninstall the BuddySuite tools.
''')

br.flags(parser, None, br.bsi_flags, None, _version)

in_args = parser.parse_args()

# Don't bother continuing on Windows. It won't work.
if sys.platform.startswith("win"):
    sys.exit("Error: Windows is not currently supported. If you really want to try BuddySuite on your current system, "
             "you will need to clone the Git repository and try accessing the modules directly. Sorry...")

if "DISPLAY" not in os.environ:  # Check that the user's system is graphical
    in_args.cmd_line = True

home_dir = path.expanduser('~')
start_dir = os.getcwd()

# Prepare working space
temp_dir = TemporaryDirectory()

source_file = path.abspath(getsourcefile(lambda: 0)).split("/")
source_file = "/".join(source_file[:-1])
with zipfile.ZipFile(source_file) as ifile:
    ifile.extractall(temp_dir.name)

os.chdir(temp_dir.name)
sys.path.insert(0, "./")

# Global variables
paths_installed = ["resources", "Bio", "dendropy", "ete3"]
files_installed = ["SeqBuddy.py", "AlignBuddy.py", "DatabaseBuddy.py", "PhyloBuddy.py", "MyFuncs.py",
                   "config.ini", "buddy_resources.py", "blastn", "blastp", "blastdbcmd"]


def check_pillow():
    try:
        from PIL import Image
        return True
    except ImportError:
        return False

# Check for pillow dependency and try to get it installed if necessary
if not in_args.cmd_line:
    if check_pillow():
        pass
    else:
        print("The image processing package 'Pillow' (https://python-pillow.github.io/) is needed for the graphical "
              "installer. If you are working within a standard desktop environment you should install Pillow, but it "
              "is not necessary if working on a non-graphical platform (e.g., on a cluster).")

        prompt = input("Try to install Pillow? (y)es/(t)emporary/[(n)o]?")
        if prompt.lower().strip() in ["yes", "y"]:
            prompt = "yes"
        elif prompt.lower().strip() in ["t", "temp", "temporary"]:
            prompt = "temp"
        else:
            prompt = "no"

        while True:  # This allows the various checks to block further checks if successful, by breakout out of loop
            if prompt in ["yes", "temp"]:
                if which("conda") and prompt == "yes":
                    print("Attempting to install Pillow with 'conda'.")
                    Popen("conda install pillow", shell=True).wait()
                    if check_pillow():
                        print("Success!")
                        break
                    else:
                        print("Failed...")

                if which("pip"):
                    print("Attempting to install Pillow with 'pip'.")
                    if prompt == "temp":
                        Popen("pip install --install-option='--prefix=%s/' --ignore-installed pillow" %
                              temp_dir.name, shell=True).wait()
                        if path.isdir("./lib/python%s/site-packages/PIL" % sys.version[:3]):
                            shutil.copytree("./lib/python%s/site-packages/PIL" % sys.version[:3], "./PIL")
                    else:
                        Popen("pip install pillow", shell=True).wait()
                    if check_pillow():
                        print("Success!")
                        break
                    else:
                        print("Failed...")

                fetch = "curl" if which("curl") else "wget" if which("wget") else False
                if fetch:
                    print("Attempting to compile Pillow from source.")
                    if fetch == "curl":
                        Popen("curl -O https://pypi.python.org/packages/source/P/Pillow/Pillow-2.9.0.tar.gz",
                              shell=True).wait()
                    else:
                        Popen("wget https://pypi.python.org/packages/source/P/Pillow/Pillow-2.9.0.tar.gz",
                              shell=True).wait()

                    Popen("tar -xzf Pillow-2.9.0.tar.gz", shell=True).wait()
                    os.chdir("Pillow-2.9.0")
                    if prompt == "temp":
                        Popen("python3 setup.py build_ext --inplace", shell=True).wait()
                        sys.path.insert(0, "%s/Pillow-2.9.0" % temp_dir.name)
                    else:
                        Popen("python3 setup.py build", shell=True).wait()
                        Popen("python3 setup.py install", shell=True).wait()
                    if check_pillow():
                        print("Success!")
                        break
                    else:
                        print("Failed...")

            print("Unable to install Pillow\nDefaulting to the command line installer.")
            in_args.cmd_line = True
            break

# Check for Tkinter
if not in_args.cmd_line:
    try:
        from tkinter import *
        from tkinter import filedialog
    except ImportError:
        proceed = False
        print("Failed to build GUI because the Tkinter package was not found. This is strange, because Tkinter is"
              "part of the Python standard library...\nDefaulting to the command line installer.")
        in_args.cmd_line = True


class BuddyInstall:
    @staticmethod
    def download_blast_binaries(install_dir, current_path, _blastn=True, _blastp=True, _blastdcmd=True):
        binary_source = 'https://raw.github.com/biologyguy/BuddySuite/master/workshop/build_dir/blast_binaries/'
        bins_to_dl = []
        current_os = sys.platform
        # Determine which binaries need to be installed
        if current_os.startswith('darwin'):
            if _blastdcmd:
                bins_to_dl.append('Darwin_blastdbcmd.zip')
            if _blastn:
                bins_to_dl.append('Darwin_blastn.zip')
            if _blastp:
                bins_to_dl.append('Darwin_blastp.zip')
        elif current_os.startswith('linux'):
            if _blastdcmd:
                bins_to_dl.append('Linux_blastdbcmd.zip')
            if _blastn:
                bins_to_dl.append('Linux_blastn.zip')
            if _blastp:
                bins_to_dl.append('Linux_blastp.zip')
        elif current_os.startswith('win'):
            if _blastdcmd:
                bins_to_dl.append('Win32_blastdbcmd.zip')
            if _blastn:
                bins_to_dl.append('Win32_blastn.zip')
            if _blastp:
                bins_to_dl.append('Win32_blastp.zip')
        else:
            return False

        file_to_name = {'Darwin_blastdbcmd.zip': 'blastdbcmd', 'Darwin_blastn.zip': 'blastn',
                        'Darwin_blastp.zip': 'blastp', 'Linux_blastdbcmd.zip': 'blastdbcmd',
                        'Linux_blastn.zip': 'blastn', 'Linux_blastp.zip': 'blastp',
                        'Win32_blastdbcmd.zip': 'blastdbcmd', 'Win32_blastn.zip': 'blastn',
                        'Win32_blastp.zip': 'blastp'}

        os.makedirs("{0}/temp".format(install_dir), exist_ok=True)
        try:
            for blast_bin in bins_to_dl:
                # Download the zip files as bytes and write them to a temporary directory
                with request.urlopen('{0}{1}'.format(binary_source, blast_bin)) as reader, \
                        open("{0}/temp/{1}".format(install_dir, blast_bin), mode='wb') as writer:
                    shutil.copyfileobj(reader, writer)
                # Unzip the zip files
                zip_file = zipfile.ZipFile("{0}/temp/{1}".format(install_dir, blast_bin))
                zip_file.extractall(path=install_dir)
                os.rename('{0}/{1}'.format(install_dir, re.sub('\.zip', '', blast_bin)),
                          '{0}/{1}'.format(install_dir, file_to_name[blast_bin]))
                # Make the binaries executable
                os.chmod('{0}/{1}'.format(install_dir, file_to_name[blast_bin]), 0o755)
                print("File added: {0}/{1}".format(current_path, file_to_name[blast_bin]))
            shutil.rmtree("{0}/temp".format(install_dir))
        except error.URLError:
            return False

        return True

    @staticmethod
    def install_buddy_suite(user_system, options):
        print("Starting.")
        BuddyInstall.edit_profile()  # Adds install directory to $PATH
        buddies_to_install = options[0]
        install_directory = options[1]
        shortcuts = options[2]
        # email_address = options[3]
        # send_diagnostics = options[4]
        # user_hash = options[5]

        for buddy in buddies_to_install:  # Deletes shortcuts from buddies that aren't being installed
            if not buddies_to_install[buddy]:
                shortcuts[buddy] = []

        if path.exists("{0}/.buddysuite/__pycache__".format(home_dir)):
            rmtree("{0}/.buddysuite/__pycache__".format(home_dir))

        uninstall = True  # Detects if BuddySuite is being uninstalled
        for buddy in buddies_to_install:
            if buddies_to_install[buddy]:
                uninstall = False

        for loc in paths_installed:  # Remove all previously installed subdirectories
            if path.exists("{0}/.buddysuite/{1}".format(home_dir, loc)):
                rmtree("{0}/.buddysuite/{1}".format(home_dir, loc))
        for loc in files_installed:  # Remove all previously installed paths
            if path.exists("{0}/.buddysuite/{1}".format(home_dir, loc)):
                os.remove("{0}/.buddysuite/{1}".format(home_dir, loc))

        # Removes the install directory if it already exists and there are no user files present
        if path.exists("{0}/.buddysuite".format(home_dir)) and len(os.listdir("{0}/.buddysuite".format(home_dir))) == 0:
            rmtree(path.realpath("{0}/.buddysuite".format(home_dir)))
            os.remove("{0}/.buddysuite".format(home_dir))

        buddy_resources_path = "./buddy_resources.py"
        myfuncs_path = "./MyFuncs.py"
        biopython_path = "./Bio"
        ete3_path = "./ete3"
        dendropy_path = "./dendropy"
        if not uninstall:  # If not uninstalling
            print("Install path: " + install_directory)
            os.makedirs(install_directory, exist_ok=True)
            print("Directory added: {0}".format(install_directory))
            if not path.exists("{0}/.buddysuite".format(home_dir)):
                os.symlink(install_directory, "{0}/.buddysuite".format(home_dir))
                print("Shortcut added: {0} ==> {1}/.buddysuite".format(install_directory, home_dir))

            if user_system in ['Darwin', 'Linux', 'Unix']:  # We don't support windows
                user_system = 'Linux' if user_system == "Unix" else user_system
                user_system = 'Win32' if user_system == "Windows" else user_system  # Not sure why this is here

                # Install all dependencies
                shutil.copy(buddy_resources_path, "{0}/buddy_resources.py".format(install_directory))
                print("File added: {0}/buddy_resources.py".format(install_directory))
                shutil.copy(myfuncs_path, "{0}/MyFuncs.py".format(install_directory))
                print("File added: {0}/MyFuncs.py".format(install_directory))
                copytree(biopython_path, "{0}/Bio".format(install_directory))
                print("Directory added: {0}/Bio".format(install_directory))
                copytree(ete3_path, "{0}/ete3".format(install_directory))
                print("Directory added: {0}/ete3".format(install_directory))
                copytree(dendropy_path, "{0}/dendropy".format(install_directory))
                print("Directory added: {0}/dendropy".format(install_directory))

                binaries = ["blastn", "blastp", "blastdbcmd"] if user_system != "Win32" else \
                    ["blastn.exe", "blastp.exe", "blastdbcmd.exe"]

                for binary in binaries:
                    if not which(binary):
                        if binary.startswith('blastn'):
                            if BuddyInstall.download_blast_binaries(install_directory, start_dir, _blastp=False, _blastdcmd=False):
                                print("{0} ==> Installed".format(binary))
                            else:
                                print("Failed to install {0}.".format(binary))
                        elif binary.startswith('blastp'):
                            if BuddyInstall.download_blast_binaries(install_directory, start_dir, _blastn=False, _blastdcmd=False):
                                print("{0} ==> Installed".format(binary))
                            else:
                                print("Failed to install {0}.".format(binary))
                        else:
                            if BuddyInstall.download_blast_binaries(install_directory, start_dir, _blastn=False, _blastp=False):
                                print("{0} ==> Installed".format(binary))
                            else:
                                print("Failed to install {0}.".format(binary))
                    else:
                        print("{0} ==> Local version detected".format(binary))

                for buddy in buddies_to_install:
                    if buddies_to_install[buddy]:
                        shutil.copy("./{0}.py".format(buddy), "{0}/{1}.py".format(install_directory, buddy))
                        os.chmod("{0}/{1}.py".format(install_directory, buddy), 0o755)
                        print("File added: {0}/{1}.py".format(install_directory, buddy))
                        for shortcut in shortcuts[buddy]:
                            if which(shortcut) is None:
                                try:
                                    os.symlink("{0}/{1}.py".format(install_directory, buddy),
                                               "{0}/{1}".format(install_directory, shortcut))
                                    print("{0}/{1}".format(install_directory, shortcut))
                                    print("Shortcut added: {0} ==> {1}".format(buddy, shortcut))
                                except FileExistsError:
                                    print("{0}/{1}".format(install_directory, shortcut))
                                    print("Shortcut retained: {0} ==> {1}".format(buddy, shortcut))

            elif user_system == 'Windows':
                print("Windows not supported at the moment.")
                return

            BuddyInstall.make_config_file(options)
        print("Finished. Please restart terminal for changes to take effect.")

    @staticmethod
    def make_config_file(options):
        print("Making config file.")
        writer = ConfigParser()
        writer.add_section('selected')
        writer.add_section('Install_path')
        writer.add_section('shortcuts')
        writer.add_section('other')
        writer['DEFAULT'] = {'selected': {'SeqBuddy': True, 'AlignBuddy': True, 'PhyloBuddy': True, 'DatabaseBuddy': True},
                             'Install_path': {'path': '{0}/BuddySuite'.format(home_dir)},
                             'shortcuts': {'SeqBuddy': 'sb\tseqbuddy', 'AlignBuddy': 'alb\talignbuddy',
                                           'PhyloBuddy': 'pb\tphylobuddy', 'DatabaseBuddy': 'db\tDatabaseBuddy'},
                             'other': {'email': '', 'diagnostics': True, 'user_hash': ''}}

        for buddy in options[0]:
            if options[0][buddy]:
                writer['selected'][buddy] = 'True'
            else:
                writer['selected'][buddy] = 'False'

        writer["Install_path"]['path'] = options[1]

        writer['other']['email'] = options[3]
        writer['other']['diagnostics'] = options[4]
        writer['other']['user_hash'] = options[5]

        for buddy in options[2]:
            sc = ''
            for shortcut in options[2][buddy]:
                sc += shortcut + "\t"
            writer['shortcuts'][buddy] = sc if sc != '' else 'None'

        with open("{0}/config.ini".format(options[1]), 'w') as configfile:
            writer.write(configfile)
        print("Config file written to {0}/config.ini".format(options[1]))

    @staticmethod
    def read_config_file():
        if path.exists("{0}/.buddysuite/config.ini".format(home_dir)):
            reader = ConfigParser()
            reader.read("{0}/.buddysuite/config.ini".format(home_dir))

            options = [{"SeqBuddy": False, "AlignBuddy": False, "PhyloBuddy": False, "DatabaseBuddy": False},
                       reader.get('Install_path', 'path'), {}, reader.get('other', 'email'),
                       reader.get('other', 'diagnostics'), reader.get('other', 'user_hash')]

            for buddy in options[0]:
                if reader['selected'][buddy] == 'True':
                    options[0][buddy] = True
                if reader['shortcuts'][buddy] != "None":
                    sc = reader['shortcuts'][buddy].split("\t")
                    options[2][buddy] = sc
                else:
                    options[2][buddy] = []

            options[4] = True if options[4] == 'True' else False

            return options

        else:
            return None

    @staticmethod
    def edit_profile():
        regex = 'export PATH=.PATH:%s/.buddysuite' % home_dir

        if system() == 'Darwin':
            if not path.exists("{0}/.profile".format(home_dir)):
                make_file = open("{0}/.profile".format(home_dir), 'w')
                make_file.close()
            with open("{0}/.profile".format(home_dir)) as file:
                contents = file.read()
                if re.search(regex, contents) is None:
                    file.close()
                    with open("{0}/.profile".format(home_dir), 'a') as file_write:
                        file_write.write("\n# added by BuddySuite installer\n")
                        file_write.write('export PATH=$PATH:{0}/.buddysuite\n'.format(home_dir))

        if system() == 'Linux':
            if not path.exists("{0}/.bashrc".format(home_dir)):
                make_file = open("{0}/.bashrc".format(home_dir), 'w')
                make_file.close()
            with open("{0}/.bashrc".format(home_dir)) as file:
                contents = file.read()
                if re.search(regex, contents) is None:
                    file.close()
                    with open("{0}/.bashrc".format(home_dir), 'a') as file_write:
                        file_write.write("\n# added by BuddySuite installer\n")
                        file_write.write('export PATH=$PATH:{0}/.buddysuite\n'.format(home_dir))


# installer
def cmd_install():
    def ask(input_prompt):
        _response = input(input_prompt)
        while True:
            if _response.lower() in ["yes", "y", '']:
                return True
            elif _response.lower() in ["no", "n", "abort"]:
                return False
            else:
                print("Response not understood. Valid options are 'yes' and 'no'.")
                _response = input(input_prompt)

    config = BuddyInstall.read_config_file()
    old_install_dir = None
    already_installed = None
    old_shortcuts = None
    install_dir = ""
    email_address = ''
    send_diagnostics = False
    user_hash = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])

    if config is not None:
        already_installed = copy.deepcopy(config[0])
        old_install_dir = copy.deepcopy(config[1])
        install_dir = old_install_dir
        old_shortcuts = copy.deepcopy(config[2])
        email_address = config[3]
        send_diagnostics = config[4]
        user_hash = config[5]

    buddies_to_install = {"SeqBuddy": False, "PhyloBuddy": False, "AlignBuddy": False, "DatabaseBuddy": False}
    shortcuts = {"SeqBuddy": [], "PhyloBuddy": [], "AlignBuddy": [], "DatabaseBuddy": []}

    def print_shortcuts():
        print("\n\tSelected Shortcuts: ")
        for buddy_tool in shortcuts:
            if buddies_to_install[buddy_tool]:
                print("\t{0}:\t\t{1}".format(buddy_tool, str(shortcuts[buddy_tool]).strip('[]')))

    def add_shortcut(in_buddy, in_shortcut):
        if buddies_to_install[in_buddy] is False:
            print("Warning: %s is not being installed. Shortcut '%s' could not be added." % (in_buddy, in_shortcut))
            return
        if which(in_shortcut) is None:
            for buddy_tool in shortcuts:
                for sc in shortcuts[buddy_tool]:
                    if in_shortcut == sc:
                        print("Warning: '{0}' is already a shortcut!".format(in_shortcut))
                        return
            shortcuts[in_buddy].append(in_shortcut)
            print("Shortcut added: {0}\t==>\t'{1}'".format(in_buddy, in_shortcut))
        else:
            print("Warning: '{0}' is already in use by another program.".format(in_shortcut))

    def remove_shortcut(in_shortcut):
        for buddy_tool in shortcuts:
            if in_shortcut in shortcuts[buddy_tool]:
                shortcuts[buddy_tool].remove(in_shortcut)
                return True
        return False

    print("\nBuddySuite command line installer launched.")
    if not ask("Would you like to proceed? ('[yes]/no') "):
        print("Installation aborted.")
        exit()

    print("\nBefore continuing, please review our license at: \nhttp://www.gnu.org/licenses/gpl-3.0.en.html")
    if not ask("Do you accept these terms? ('[yes]/no') "):
        print("Installation aborted.")
        exit()

    if config is not None:
        print("\nWe have detected a previous installation. Some settings will be imported.")

    for buddy in buddies_to_install:
        operation = ' install' if already_installed is None or not already_installed[buddy] else ' keep'
        optional = '' if already_installed is None or not already_installed[buddy] else ' installed'
        if ask("Would you like to{0} {1}{2}? ('[yes]/no') ".format(operation, buddy, optional)):
            buddies_to_install[buddy] = True
        else:
            buddies_to_install[buddy] = False

    all_false = True
    for buddy in buddies_to_install:
        if buddies_to_install[buddy]:
            all_false = False

    if not all_false:
        if old_shortcuts is not None:
            print("\nThe following shortcuts already exist on your system:")
            for buddy in old_shortcuts:
                if buddies_to_install[buddy]:
                    print("\t{0}:\t{1}".format(buddy, str(old_shortcuts[buddy]).strip('[]')))

            if ask("\nWould you like to keep these old shortcuts? ('[yes]/no') "):
                shortcuts = old_shortcuts

        default_shortcuts = {"SeqBuddy": ['sb', 'seqbuddy'], "PhyloBuddy": ['pb', 'phylobuddy'],
                             "AlignBuddy": ['alb', 'alignbuddy'], "DatabaseBuddy": ['db', 'dbbuddy']}

        for buddy in default_shortcuts:
            for shortcut in default_shortcuts[buddy]:
                if which(shortcut) is not None:
                    default_shortcuts[buddy].remove(shortcut)

        for _buddy in shortcuts:
            for _shortcut in shortcuts[_buddy]:
                for buddy in default_shortcuts:
                    for shortcut in default_shortcuts[buddy]:
                        if shortcut == _shortcut:
                            default_shortcuts[buddy].remove(shortcut)
        empty = True
        for buddy in default_shortcuts:
            if buddies_to_install[buddy] and len(default_shortcuts[buddy]) > 0:
                empty = False
        if not empty:
            print("\tDefault Shortcuts:")
            for buddy in default_shortcuts:
                if buddies_to_install[buddy]:
                    print("\t{0}:\t{1}".format(buddy, str(default_shortcuts[buddy]).strip('[]')))

            if ask("\nWould you like to add the default shortcuts? ('[yes]/no') "):
                for buddy in default_shortcuts:
                    if buddies_to_install[buddy]:
                        for shortcut in default_shortcuts[buddy]:
                            add_shortcut(buddy, shortcut)

        print_shortcuts()

        if ask("\nWould you like to add or remove shortcuts? ('[yes]/no') "):
            _prompt = True
            while _prompt:
                add_or_remove = input("Would you like to add or remove? ('add/remove/[cancel]') ")
                while True:
                    if add_or_remove.lower() in ['cancel', '']:
                        print("Cancelled.")
                        add_or_remove = None
                        break
                    elif add_or_remove.lower() in ['add', 'a']:
                        add_or_remove = 'add'
                        break
                    elif add_or_remove.lower() in ['remove', 'r', 'rm']:
                        add_or_remove = 'remove'
                        break
                    else:
                        print("Response not understood. Try again.")
                        add_or_remove = input("Would you like to add or remove? ('add/remove/cancel') ")

                if add_or_remove == 'add':
                    response = input("What would you like to add? (shortcut name) ")
                    response = re.sub("[^a-zA-Z0-9]", '', response)
                    buddy_string = ''
                    for buddy in buddies_to_install:
                        if buddies_to_install[buddy]:
                            buddy_string += "{0}/".format(buddy)
                    buddy_string += '[cancel]'
                    buddy = input("What are you linking to? ('{0}') ".format(buddy_string))
                    while True:
                        if buddy.lower() in ['cancel', '']:
                            print("Shortcut not added: {0}".format(response))
                            break
                        elif (len(buddy) >= 2 and buddy in 'seq') or buddy.lower() == 'seqbuddy':
                            add_shortcut("SeqBuddy", response)
                            break
                        elif (len(buddy) >= 2 and buddy in 'phylo') or buddy.lower() == 'phylobuddy':
                            add_shortcut("PhyloBuddy", response)
                            break
                        elif (len(buddy) >= 2 and buddy in 'align') or buddy.lower() == 'alignbuddy':
                            add_shortcut("AlignBuddy", response)
                            break
                        elif (len(buddy) >= 2 and buddy in 'data') or buddy.lower() == 'databasebuddy':
                            add_shortcut("DatabaseBuddy", response)
                            break
                        else:
                            print("Response not understood. Try again.")
                            buddy = input("What are you linking to? "
                                          "('{0}') ".format(buddy_string))
                elif add_or_remove == 'remove':
                    response = input("What would you like to remove? (shortcut name) ")
                    if remove_shortcut(response):
                        print("Shortcut removed: '{0}'".format(response))
                    else:
                        print("Warning: Shortcut already does not exist.")

                print_shortcuts()
                _prompt = ask("Would you like to add/remove another shortcut? ('[yes]/no') ")

        if old_install_dir is None:
            install_dir = input("\nBuddySuite will be installed in your home directory, or specify a different "
                                "location ('/BuddySuite' will be appended): ")
            install_dir = install_dir.rstrip("/")
            while True:
                if install_dir.lower() == 'abort':
                    print('Installation aborted.')
                    exit()
                try:
                    install_dir = re.sub('~', home_dir, install_dir)
                    if install_dir == "":
                        install_dir = home_dir
                    else:
                        install_dir = os.path.abspath(install_dir)
                    install_dir = install_dir.rstrip("/")
                    os.makedirs("{0}/BuddySuite".format(install_dir), exist_ok=True)
                    install_dir = "{0}/BuddySuite".format(install_dir)
                    break
                except PermissionError:
                    print('Insufficient privileges to write here.')

                install_dir = input("Please specify a different installation directory (or say 'abort') ")
        else:
            install_dir = old_install_dir

    print("\nProviding a valid email address is recommended if using DatabaseBuddy to access NCBI, as they will attempt "
          "to contact you before blocking your IP if you are not adhering to their usage limitations. "
          "See our privacy statement for further details.\n"
          "https://github.com/biologyguy/BuddySuite/blob/master/privacy\n")

    if email_address != '':
        if not ask("Your email address is currently set to \033[1m%s\033[m, would you like to keep it the same? [y]/n"
                   % email_address):
            email_address = ''

    if email_address == '':
        email_address = input("\nPlease provide your email address (optional): ")

    if ask("\nWould you like to join our Software Improvement Program (SIP)? Anonymized usage statistics and crash "
           "reports will be automatically transmitted to the BuddySuite developers. [y]/n"):
        send_diagnostics = True

    print("\nPlease verify your settings.\n")
    for buddy in buddies_to_install:
        if already_installed is not None and already_installed[buddy]:
            if buddies_to_install[buddy]:
                print("\t{0}:\t\tModify".format(buddy))
            else:
                print("\t{0}:\t\tUninstall".format(buddy))
        else:
            if buddies_to_install[buddy]:
                print("\t{0}:\tInstall".format(buddy))
            else:
                print("\t{0}:\tSkip".format(buddy))
    if not all_false:
        print_shortcuts()
        print("\n\tInstallation directory:\t{0}\n".format(install_dir))

    print("\tEmail address:\t\t\t%s" % email_address)
    print("\tSIP participation:\t%s\n" % send_diagnostics)

    if ask("Are these settings okay? ('[yes]/abort') "):
        if all_false and config is not None:
            os.remove("{0}/.buddysuite/config.ini".format(home_dir))
        if config is not None:  # Uninstall removed shortcuts
            for buddy in old_shortcuts:
                for shortcut in old_shortcuts[buddy]:
                    if os.path.exists("%s/.buddysuite/%s" % (home_dir, shortcut)) and shortcut not in shortcuts[buddy]:
                        os.remove("{0}/.buddysuite/{1}".format(home_dir, shortcut))
        for buddy in buddies_to_install:
            if not buddies_to_install[buddy]:
                shortcuts[buddy] = []
        BuddyInstall.install_buddy_suite(system(), [buddies_to_install, install_dir, shortcuts,
                                                    email_address, str(send_diagnostics), user_hash])
        exit()
    else:
        print("Installation aborted.")
        exit()

if in_args.cmd_line:
    cmd_install()
else:
    root = Tk()
    sw = root.winfo_screenwidth()
    sh = root.winfo_screenheight()
    scale_factor = sh / 1440
    # print(scale_factor)
    sys.path.insert(0, "./")
    root.title("BuddySuite Installer")
    root.config(background='white')
    root.tk_setPalette(background='white')


class Installer(Frame):
    from PIL import Image
    container = []
    buddies = collections.OrderedDict()
    buddy_names = ["SeqBuddy", "PhyloBuddy", "AlignBuddy", "DatabaseBuddy"]
    for buddy in buddy_names:
        buddies[buddy] = True

    uninstall = False
    image_paths = ["{0}/BuddySuite-logo.gif".format(temp_dir.name), "{0}/InstallDirectory.gif".format(temp_dir.name),
                   "{0}/ConsoleShortcuts.gif".format(temp_dir.name), "{0}/ConfirmSelection.gif".format(temp_dir.name),
                   "{0}/AlignBuddy-logo.gif".format(temp_dir.name), "{0}/SeqBuddy-logo.gif".format(temp_dir.name),
                   "{0}/PhyloBuddy-logo.gif".format(temp_dir.name), "{0}/DBBuddy-logo.gif".format(temp_dir.name),
                   "{0}/EmailAddress.gif".format(temp_dir.name), "{0}/SendDiagnostics.gif".format(temp_dir.name)]

    for _path in image_paths:
        image = Image.open(_path)
        new_width = image.size[0] * scale_factor
        new_height = image.size[1] * scale_factor
        new_size = new_width, new_height
        image.thumbnail(new_size, Image.ANTIALIAS)
        image.save(_path, "gif")

    bs_logo = PhotoImage(file="{0}/BuddySuite-logo.gif".format(temp_dir.name))
    id_logo = PhotoImage(file="{0}/InstallDirectory.gif".format(temp_dir.name))
    sc_logo = PhotoImage(file="{0}/ConsoleShortcuts.gif".format(temp_dir.name))
    cs_logo = PhotoImage(file="{0}/ConfirmSelection.gif".format(temp_dir.name))
    ea_logo = PhotoImage(file="{0}/EmailAddress.gif".format(temp_dir.name))
    sd_logo = PhotoImage(file="{0}/SendDiagnostics.gif".format(temp_dir.name))

    alb_logo = PhotoImage(file="{0}/AlignBuddy-logo.gif".format(temp_dir.name))
    sb_logo = PhotoImage(file="{0}/SeqBuddy-logo.gif".format(temp_dir.name))
    pb_logo = PhotoImage(file="{0}/PhyloBuddy-logo.gif".format(temp_dir.name))
    db_logo = PhotoImage(file="{0}/DBBuddy-logo.gif".format(temp_dir.name))

    suite_logos = [sb_logo, pb_logo, alb_logo, db_logo]

    email_address = ''
    send_diagnostics = True
    user_hash = "".join([choice(string.ascii_letters + string.digits) for _ in range(10)])

    install_dir = "{0}/BuddySuite".format(home_dir)
    default_dir = "{0}/BuddySuite".format(home_dir)
    default = True
    shortcuts = {"SeqBuddy": ['sb', 'seqbuddy'], "PhyloBuddy": ['pb', 'phylobuddy'],
                 "AlignBuddy": ['alb', 'alignbuddy'], "DatabaseBuddy": ['db', 'dbbuddy']}
    for buddy in shortcuts:
        for shortcut in shortcuts[buddy]:
            if which(shortcut) is not None:
                shortcuts[buddy].remove(shortcut)
    user_system = system()
    user_os = platform()

    modifying = False
    config = None

    if BuddyInstall.read_config_file():
        modifying = True
        config = BuddyInstall.read_config_file()
        buddies = copy.deepcopy(config[0])
        install_dir = copy.deepcopy(config[1])
        shortcuts = copy.deepcopy(config[2])
        email_address = config[3]
        send_diagnostics = config[4]
        user_hash = config[5]

        # if new install of given tool, re-add default shortcuts]
        if which("sb") is None and not buddies["SeqBuddy"]:
            shortcuts["SeqBuddy"].append("sb")
        if which("seqbuddy") is None and not buddies["SeqBuddy"]:
            shortcuts["SeqBuddy"].append("seqbuddy")
        if which("alb") is None and not buddies["AlignBuddy"]:
            shortcuts["AlignBuddy"].append("alb")
        if which("alignbuddy") is None and not buddies["AlignBuddy"]:
            shortcuts["AlignBuddy"].append("alignbuddy")
        if which("pb") is None and not buddies["PhyloBuddy"]:
            shortcuts["PhyloBuddy"].append("pb")
        if which("phylobuddy") is None and not buddies["PhyloBuddy"]:
            shortcuts["PhyloBuddy"].append("phylobuddy")
        if which("db") is None and not buddies["DatabaseBuddy"]:
            shortcuts["DatabaseBuddy"].append("db")
        if which("databasebuddy") is None and not buddies["DatabaseBuddy"]:
            shortcuts["DatabaseBuddy"].append("dbbuddy")

    original_shortcuts = copy.deepcopy(shortcuts)
    conflict = False

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.welcome()
        self.pack()

    def uninstall_all(self):
        self.uninstall = True
        for buddy in self.buddies:
            self.buddies[buddy] = False
        self.confirmation()

    def welcome(self, *args):
        self.clear_container()
        title_frame = Frame()
        welcome_label = Label(title_frame, image=self.bs_logo)
        welcome_label.pack(side=TOP)
        version_label = Label(title_frame, text="Version {0}".format(_version.short()))
        version_label.pack(side=RIGHT)
        title_frame.pack(pady=sh / 10)
        self.container.append(title_frame)
        button_container = Frame()
        next_button = Button(button_container, width=20, pady=int(20 * scale_factor),
                             text="Install", command=self.license)
        uninstall_button = Button(button_container, width=20, pady=int(20 * scale_factor),
                                  text="Uninstall", command=self.uninstall_all)
        if self.modifying:
            next_button.config(text="Modify Installation")
            uninstall_button.pack(side=BOTTOM)
        next_button.pack(side=TOP)
        button_container.pack(side=BOTTOM, pady=40 * scale_factor)
        self.container.append(button_container)

    def license(self):
        self.uninstall = False
        self.clear_container()
        frame = Frame(pady=75 * scale_factor)
        scrollbar = Scrollbar(master=frame)
        license_file = open("{0}/LICENSE".format(temp_dir.name))
        license_box = Text(master=frame, wrap=WORD, yscrollcommand=scrollbar.set, width=int(75 * scale_factor),
                           height=int(25 * scale_factor))
        license_box.insert(END, license_file.read())
        license_box.config(state=DISABLED)
        scrollbar.config(command=license_box.yview)

        frame.pack(side=TOP)
        license_box.pack(side=LEFT)
        scrollbar.pack(side=RIGHT, fill=Y)
        button_frame = Frame()
        next_button = Button(button_frame, padx=50, pady=20 * scale_factor, text="I agree", command=self.next_tool)
        next_button.pack(side=TOP)
        back_button = Label(button_frame, padx=50, pady=20 * scale_factor, text="Cancel", fg='blue',
                            font=('Helvetica', 12, 'underline'))
        back_button.bind("<Button-1>", self.welcome)
        back_button.pack(side=BOTTOM)
        button_frame.pack(side=BOTTOM, pady=20 * scale_factor)
        self.container.append(frame)
        self.container.append(button_frame)

    def next_tool(self, num=0, entry=None):
        if entry:
            self.install_dir = entry.get()

        self.clear_container()

        all_false = True
        for buddy in self.buddies:
            if self.buddies[buddy]:
                all_false = False
        if num < 0:
            self.license()
            return
        elif num > 3:
            if all_false:
                if self.modifying:
                    self.confirmation()
                else:
                    self.none_selected_page()
            else:
                self.install_location()
            return

        logo_label = Label(image=self.suite_logos[num], pady=20 * scale_factor)
        logo_label.pack(side=TOP)
        self.container.append(logo_label)
        mega_frame = Frame(pady=20 * scale_factor)
        self.container.append(mega_frame)
        frame = Frame(mega_frame, padx=50 * scale_factor, pady=10 * scale_factor)
        scrollbar = Scrollbar(master=frame)
        description_box = Text(master=frame, wrap=WORD, yscrollcommand=scrollbar.set, width=int(75 * scale_factor))
        if num == 0:
            description_box.insert(END, seqbuddy_blurb)
        elif num == 1:
            description_box.insert(END, phylobuddy_blurb)
        elif num == 2:
            description_box.insert(END, alignbuddy_blurb)
        elif num == 3:
            description_box.insert(END, databasebuddy_blurb)
        else:
            description = open("{0}/LICENSE".format(temp_dir.name)).read()
            description_box.insert(END, description)

        description_box.config(state=DISABLED)
        scrollbar.config(command=description_box.yview)
        description_box.pack(side=LEFT)
        scrollbar.pack(side=RIGHT, fill=Y)
        button_frame = Frame(mega_frame)
        next_func = partial(self.next_tool, num + 1)
        next_button = Button(button_frame, padx=50 * scale_factor,
                             pady=20 * scale_factor, text="Next", command=next_func)
        self.container.append(next_button)
        next_button.pack(side=RIGHT)
        back_func = partial(self.next_tool, num - 1)
        back_button = Button(button_frame, padx=50 * scale_factor,
                             pady=20 * scale_factor, text="Back", command=back_func)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM)
        _func = partial(self.toggle_tool, self.buddy_names[num])
        if not self.modifying or self.config[0][self.buddy_names[num]] is False:
            tool_button = Checkbutton(mega_frame, text="Install {0}".format(self.buddy_names[num]), command=_func,
                                      pady=20 * scale_factor)
            if self.buddies[self.buddy_names[num]]:
                tool_button.select()
            else:
                tool_button.deselect()
            tool_button.pack(side=BOTTOM)
        else:
            radio_frame = Frame(mega_frame)
            var = IntVar()
            if self.buddies[self.buddy_names[num]]:
                var.set(1)
            else:
                var.set(0)
            radiobutton_func = partial(self.toggle_tool, self.buddy_names[num], var)
            update = Radiobutton(radio_frame, text="Update/Repair", value=1, variable=var, command=radiobutton_func,
                                 pady=20 * scale_factor)
            uninstall = Radiobutton(radio_frame, text="Uninstall", value=0, variable=var, command=radiobutton_func,
                                    pady=20 * scale_factor)
            update.pack(side=LEFT)
            uninstall.pack(side=RIGHT)
            radio_frame.pack(side=BOTTOM)
        frame.pack(side=TOP)
        mega_frame.pack(side=TOP)

    def toggle_tool(self, name, var=None):
        if var is not None:
            if var.get() == 1:
                self.buddies[name] = True
            else:
                self.buddies[name] = False
        else:
            self.buddies[name] = False if self.buddies[name] else True

    def none_selected_page(self):
        self.clear_container()
        labelframe = Frame(root)
        label1 = Label(labelframe, text="You're currently not installing anything!", pady=5 * scale_factor,
                       font=('Helvetica', int(16 * scale_factor)))
        label2 = Label(labelframe, text="Please go back and select at least one tool.", pady=5 * scale_factor,
                       font=('Helvetica', int(16 * scale_factor)))
        label1.pack(side=TOP)
        label2.pack(side=BOTTOM)
        labelframe.pack(side=TOP, pady=120 * scale_factor)
        back_func = partial(self.next_tool, 3)
        back_button = Button(root, padx=50 * scale_factor, pady=20 * scale_factor, text="Back", command=back_func)
        back_button.pack(side=BOTTOM, pady=40 * scale_factor)
        self.container.append(labelframe)
        self.container.append(back_button)

    def install_location(self, email=None):
        if email is not None:
            self.email_address = email.get()
        self.clear_container()
        logo_label = Label(image=self.id_logo, pady=20 * scale_factor)
        logo_label.pack(side=TOP)
        self.container.append(logo_label)
        frame = Frame(root, padx=50 * scale_factor, pady=50 * scale_factor)
        label = Label(frame, text="Directory:")
        label.pack(side=TOP, anchor=NW)
        directory_frame = Frame(frame)
        directory_text = Entry(directory_frame)
        directory_text.insert(END, self.install_dir)
        browse_func = partial(self.choose_directory, directory_text)
        browse_button = Button(directory_frame, text="Browse", command=browse_func)
        browse_button.pack(side=RIGHT)
        directory_text.pack(side=TOP, anchor=NW, fill=X)
        directory_frame.pack(side=TOP, anchor=NW, fill=X)
        frame.pack(side=TOP, padx=10 * scale_factor, pady=10 * scale_factor, fill=BOTH)

        toggle_func = partial(self.default_directory, directory_text, browse_button)
        toggle_default = Checkbutton(frame, text="Default directory", pady=10 * scale_factor, command=toggle_func)
        toggle_default.pack(side=LEFT)

        if self.default:
            browse_button.config(state=DISABLED)
            directory_text.config(state=DISABLED)
            toggle_default.select()
        else:
            toggle_default.deselect()

        self.container.append(frame)
        lower_box = Frame()
        if self.modifying:
            browse_button.config(state=DISABLED)
            directory_text.config(state=DISABLED)
            toggle_default.config(state=DISABLED)
            warning = Label(lower_box, text="Previous install detected. Uninstall first to change install directory.")
            warning.pack(side=TOP, pady=50 * scale_factor)
        button_frame = Frame(lower_box)
        next_func = partial(self.email_address_page, directory_text)
        next_button = Button(button_frame, padx=50 * scale_factor,
                             pady=20 * scale_factor, text="Next", command=next_func)
        next_button.pack(side=RIGHT)
        back_func = partial(self.next_tool, 3, directory_text)
        back_button = Button(button_frame, padx=50 * scale_factor,
                             pady=20 * scale_factor, text="Back", command=back_func)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM, pady=20 * scale_factor)
        lower_box.pack(side=BOTTOM)
        self.container.append(lower_box)

    def email_address_page(self, in_dir=None):
        if in_dir is not None:
            self.install_dir = in_dir.get()
        self.clear_container()

        logo_label = Label(image=self.ea_logo, pady=20 * scale_factor)
        logo_label.pack(side=TOP)
        self.container.append(logo_label)
        info_frame = LabelFrame(bd=0, padx=10 * scale_factor, pady=10 * scale_factor)
        self.container.append(info_frame)
        text_lable = Label(info_frame, justify="left", wraplength=700 * scale_factor,
                           text="Providing a valid email address is recommended if using DatabaseBuddy to access NCBI, "
                                "as they will attempt to contact you before blocking your IP if you are not adhering "
                                "to their usage limitations. See our privacy statement for further details.\n"
                                "https://github.com/biologyguy/BuddySuite/blob/master/privacy")
        text_lable.grid(row=0, sticky=NW)
        info_frame.pack(side=TOP, anchor=NW, padx=50 * scale_factor, pady=50 * scale_factor, fill=BOTH)

        frame = Frame(root, padx=50 * scale_factor, pady=50 * scale_factor)
        label = Label(frame, text="Please provide your email address (optional).")
        label.pack(side=TOP, anchor=NW)
        email_text = Entry(frame)
        email_text.insert(END, self.email_address)
        email_text.pack(side=TOP, anchor=NW, fill=X)

        frame.pack(side=TOP, padx=10 * scale_factor, pady=10 * scale_factor, fill=BOTH)
        self.container.append(frame)

        lower_box = Frame()

        button_frame = Frame(lower_box)
        next_func = partial(self.install_shortcuts, email_text)
        next_button = Button(button_frame, padx=50 * scale_factor,
                             pady=20 * scale_factor, text="Next", command=next_func)
        next_button.pack(side=RIGHT)
        back_func = partial(self.install_location, email_text)
        back_button = Button(button_frame, padx=50 * scale_factor,
                             pady=20 * scale_factor, text="Back", command=back_func)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM, pady=20 * scale_factor)
        lower_box.pack(side=BOTTOM)
        self.container.append(lower_box)

    def install_shortcuts(self, email=None):
        if email is not None:
            self.email_address = email.get()
        self.clear_container()
        logo_label = Label(image=self.sc_logo)
        logo_label.pack(side=TOP)
        self.container.append(logo_label)

        button_frame = Frame()
        next_button = Button(button_frame, padx=50 * scale_factor, pady=20 * scale_factor, text="Next",
                             command=self.diagnostic_prompt)
        next_button.pack(side=RIGHT)
        back_button = Button(button_frame, padx=50 * scale_factor, pady=20 * scale_factor, text="Back",
                             command=self.email_address_page)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM, pady=40 * scale_factor)
        self.container.append(button_frame)

        frame = Frame()
        self.container.append(frame)
        scrollbox_frame = Frame(frame)
        scrollbar = Scrollbar(master=scrollbox_frame)
        shortcut_box = Listbox(master=scrollbox_frame, yscrollcommand=scrollbar.set, bd=2, relief=SUNKEN,
                               height=int(12 * scale_factor))
        for buddy in self.shortcuts:
            if self.buddies[buddy]:
                for shortcut in self.shortcuts[buddy]:
                    shortcut_box.insert(END, "{0} ==> {1}".format(buddy, shortcut))
        scrollbar.config(command=shortcut_box.yview)
        shortcut_box.pack(side=LEFT, fill=BOTH, expand=1)
        scrollbar.pack(side=RIGHT, fill=Y)
        scrollbox_frame.pack(side=TOP, fill=BOTH, expand=1)
        space = Frame()
        space.pack(pady=25 * scale_factor)

        self.container.append(space)

        debug_frame = Frame(frame)
        debug = Label(debug_frame)
        debug.pack(side=BOTTOM, anchor=NW)
        selection_frame = Frame(debug_frame)
        entry_frame = Frame(selection_frame)
        buds = [x for x in self.buddies if self.buddies[x] is True]
        curr_buddy = StringVar()
        for buddy in self.buddies:
            if self.buddies[buddy]:
                curr_buddy.set(buddy)
                break
        dropdown = OptionMenu(selection_frame, curr_buddy, *buds)
        dropdown.config(width=13)
        dropdown.pack(side=RIGHT)
        selection_frame.pack(fill=X, expand=1, side=TOP)
        debug_frame.pack(fill=X, expand=1, side=BOTTOM, anchor=NW)
        self.container.append(entry_frame)
        entry_button_frame = Frame(entry_frame)
        shortcut_entry = Entry(entry_frame, width=int(20 * scale_factor))
        add_func = partial(self.add_shortcut, curr_buddy, shortcut_box, shortcut_entry, debug)
        add_button = Button(entry_button_frame, text="Add", command=add_func, padx=1)
        add_button.pack(side=LEFT)
        rmv_func = partial(self.remove_shortcut, curr_buddy, shortcut_box, shortcut_entry, debug)
        rmv_button = Button(entry_button_frame, text="Remove", command=rmv_func, padx=1)
        rmv_button.pack(side=RIGHT)
        shortcut_entry.pack(side=LEFT, fill=X, expand=1)
        entry_button_frame.pack(side=LEFT)
        frame.pack(padx=100 * scale_factor, expand=1, fill=BOTH, side=BOTTOM)
        entry_frame.pack(fill=X, expand=1, side=LEFT)

        click_func = partial(self.click_shortcut, shortcut_entry, curr_buddy)
        shortcut_box.bind("<Button-1>", click_func)

    @staticmethod
    def click_shortcut(entry, var, event):
        event.widget.activate(event.widget.nearest(event.y))
        entry.delete(0, END)
        text = event.widget.get(ACTIVE).split(" ==> ")
        var.set(text[0])
        entry.insert(END, text[1])

    def add_shortcut(self, buddy, listbox, entry, debug):
        debug.config(text='')
        addable = True
        text = re.sub("[^a-zA-Z0-9]", '', entry.get())
        entry.delete(0, END)
        entry.insert(END, text)
        bud = buddy.get()
        for name in self.buddy_names:
            if "{0} ==> {1}".format(name, text) in listbox.get(0, END):
                addable = False
        if addable:
            if which(text) is not None and text not in self.original_shortcuts[bud]:
                debug.config(text='Not added: Naming Conflict')
            else:
                self.shortcuts[bud].append(text)
                listbox.delete(0, END)
                for buddy in self.shortcuts:
                    if self.buddies[buddy]:
                        for shortcut in self.shortcuts[buddy]:
                            listbox.insert(END, "{0} ==> {1}".format(buddy, shortcut))
                debug.config(text='')

        else:
            debug.config(text='Not added: Already exists')
            return

    def remove_shortcut(self, buddy, listbox, entry, debug):
        debug.config(text='')
        text = re.sub("[^a-zA-Z0-9]", '', entry.get())
        entry.delete(0, END)
        entry.insert(END, text)
        lst = listbox.get(0, END)
        if "{0} ==> {1}".format(buddy.get(), text) in lst:
            index = lst.index("{0} ==> {1}".format(buddy.get(), text))
            listbox.delete(index)
            self.shortcuts[buddy.get()].remove(text)

    def diagnostic_prompt(self):
        self.clear_container()
        logo_label = Label(image=self.sd_logo, pady=20 * scale_factor)
        logo_label.pack(side=TOP)
        self.container.append(logo_label)

        frame = Frame(root, padx=50 * scale_factor, pady=75 * scale_factor)

        diagnostic_box = Checkbutton(frame,
                                     text="  Please join our Software Improvement Program.\n"
                                          "  Anonymized usage statistics and crash reports will be\n"
                                          "  automatically transmitted to the BuddySuite developers.",
                                     command=self.toggle_diagnostics, pady=50 * scale_factor, justify="left")
        if self.send_diagnostics:
            diagnostic_box.toggle()

        diagnostic_box.pack(side=BOTTOM)

        frame.pack(side=TOP, padx=10 * scale_factor, pady=10 * scale_factor, fill=BOTH)
        self.container.append(frame)

        lower_box = Frame()

        button_frame = Frame(lower_box)
        next_button = Button(button_frame, padx=50 * scale_factor,
                             pady=20 * scale_factor, text="Next", command=self.confirmation)
        next_button.pack(side=RIGHT)
        back_button = Button(button_frame, padx=50 * scale_factor,
                             pady=20 * scale_factor, text="Back", command=self.install_shortcuts)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM, pady=20 * scale_factor)
        lower_box.pack(side=BOTTOM)
        self.container.append(lower_box)

    def toggle_diagnostics(self):
        if self.send_diagnostics:
            self.send_diagnostics = False
        else:
            self.send_diagnostics = True

    def confirmation(self):
        confirmation_font = ('Courier', int(15 * scale_factor))
        self.clear_container()
        logo_label = Label(image=self.cs_logo, pady=18 * scale_factor)
        logo_label.pack(side=TOP)
        self.container.append(logo_label)
        info_frame = LabelFrame(text="Selections", bd=2, relief=SUNKEN, padx=10 * scale_factor, pady=10 * scale_factor)
        self.container.append(info_frame)
        os_label = Label(info_frame, text="{:<28}{}".format("Operating System:", self.user_os),
                         font=confirmation_font)
        dir_label = Label(info_frame, text="{:<28}{}".format("Install Directory:", self.install_dir),
                          font=confirmation_font)

        if not self.modifying or self.config[0]["SeqBuddy"] is False:
            install_state = "Install" if self.buddies["SeqBuddy"] else "Skip"
        else:
            install_state = "Modify" if self.buddies["SeqBuddy"] else "Uninstall"
        sb_label = Label(info_frame, text="{:<28}{}".format("SeqBuddy:", install_state),
                         font=confirmation_font)

        if not self.modifying or self.config[0]["PhyloBuddy"] is False:
            install_state = "Install" if self.buddies["PhyloBuddy"] else "Skip"
        else:
            install_state = "Modify" if self.buddies["PhyloBuddy"] else "Uninstall"
        pb_label = Label(info_frame, text="{:<28}{}".format("PhyloBuddy:", install_state), font=confirmation_font)

        if not self.modifying or self.config[0]["AlignBuddy"] is False:
            install_state = "Install" if self.buddies["AlignBuddy"] else "Skip"
        else:
            install_state = "Modify" if self.buddies["AlignBuddy"] else "Uninstall"
        ab_label = Label(info_frame, text="{:<28}{}".format("AlignBuddy:", install_state),
                         font=confirmation_font)

        if not self.modifying or self.config[0]["DatabaseBuddy"] is False:
            install_state = "Install" if self.buddies["DatabaseBuddy"] else "Skip"
        else:
            install_state = "Modify" if self.buddies["DatabaseBuddy"] else "Uninstall"
        db_label = Label(info_frame, text="{:<28}{}".format("DatabaseBuddy:", install_state),
                         font=confirmation_font)

        cs_sb_label = Label(info_frame, font=confirmation_font,
                            text="{:<28}{}".format("SeqBuddy Shortcut(s):", ", ".join(self.shortcuts["SeqBuddy"])))
        cs_ab_label = Label(info_frame, font=confirmation_font,
                            text="{:<28}{}".format("AlignBuddy Shortcut(s):", ", ".join(self.shortcuts["AlignBuddy"])))
        cs_pb_label = Label(info_frame, font=confirmation_font,
                            text="{:<28}{}".format("PhyloBuddy Shortcut(s):", ", ".join(self.shortcuts["PhyloBuddy"])))
        cs_db_label = Label(info_frame, font=confirmation_font,
                            text="{:<28}{}".format("DatabaseBuddy Shortcut(s):", ", ".join(self.shortcuts["DatabaseBuddy"])))

        os_label.grid(row=0, sticky=NW)
        dir_label.grid(row=1, sticky=NW)
        sb_label.grid(row=2, sticky=NW)
        pb_label.grid(row=3, sticky=NW)
        ab_label.grid(row=4, sticky=NW)
        db_label.grid(row=5, sticky=NW)

        if self.buddies["SeqBuddy"]:
            cs_sb_label.grid(row=6, sticky=NW)
        if self.buddies["AlignBuddy"]:
            cs_ab_label.grid(row=7, sticky=NW)
        if self.buddies["PhyloBuddy"]:
            cs_pb_label.grid(row=8, sticky=NW)
        if self.buddies["DatabaseBuddy"]:
            cs_db_label.grid(row=9, sticky=NW)

        if self.email_address != '' and self.buddies["DatabaseBuddy"]:
            em_label = Label(info_frame, font=confirmation_font,
                             text="Email address: {0}".format(self.email_address))
            em_label.grid(row=10, sticky=NW)

        dp_label = Label(info_frame, font=confirmation_font,
                         text="Send Diagnostics: {0}".format(self.send_diagnostics))
        dp_label.grid(row=11, sticky=NW)

        info_frame.pack(side=TOP, anchor=NW, padx=50 * scale_factor, pady=50 * scale_factor, fill=BOTH)

        all_false = True
        for buddy in self.buddies:
            if self.buddies[buddy]:
                all_false = False
        if self.config is not None and all_false:
            back_func = partial(self.next_tool, 3)
        else:
            back_func = self.diagnostic_prompt

        button_frame = Frame()
        next_button = Button(button_frame, padx=50 * scale_factor, pady=20 * scale_factor, text="Install",
                             command=self.install)
        if all_false:
            next_button.config(text="Uninstall")
        next_button.pack(side=RIGHT)
        back_button = Button(button_frame, padx=50 * scale_factor,
                             pady=20 * scale_factor, text="Back", command=back_func)
        if self.uninstall:
            back_button.config(command=self.welcome)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM, pady=40 * scale_factor)
        self.container.append(button_frame)

    def install(self):
        all_false = True

        for buddy in self.buddies:
            if self.buddies[buddy]:
                all_false = False
            else:
                self.shortcuts[buddy] = []

        if all_false and self.config is not None:
            os.remove("{0}/.buddysuite/config.ini".format(home_dir))
        if self.config is not None:  # Uninstall removed shortcuts
            for buddy in self.original_shortcuts:
                for shortcut in self.original_shortcuts[buddy]:
                    if os.path.exists("{0}/.buddysuite/{1}".format(home_dir, shortcut)) \
                            and shortcut not in self.shortcuts[buddy]:
                        os.remove("{0}/.buddysuite/{1}".format(home_dir, shortcut))
        for buddy in self.buddies:
            if not self.buddies[buddy]:
                self.shortcuts[buddy] = []
        BuddyInstall.install_buddy_suite(self.user_system, [self.buddies, self.install_dir, self.shortcuts,
                                                            self.email_address, str(self.send_diagnostics),
                                                            self.user_hash])
        exit()

    def clear_container(self):
        for item in self.container:
            item.destroy()

    @staticmethod
    def choose_directory(textbox):
        name = filedialog.askdirectory(parent=root, title="Select Directory", initialdir=textbox.get())
        if name is not "":
            textbox.delete(0, len(textbox.get()))
            textbox.insert(END, name + "/BuddySuite")

    def default_directory(self, textbox, browse):
        if self.default:
            textbox.config(state=NORMAL)
            browse.config(state=NORMAL)
            self.default = False
        else:
            self.install_dir = self.default_dir
            textbox.delete(0, len(textbox.get()))
            textbox.insert(END, self.default_dir)
            browse.config(state=DISABLED)
            textbox.config(state=DISABLED)
            self.default = True

# ################################################### Buddy Blurbs ################################################### #
seqbuddy_blurb = "Read, write, analyze, and manipulate sequence files in common formats including FASTA, " \
                 "GenBank, and EMBL. From the command line, SeqBuddy input can be file paths or piped data " \
                 "and the content is automatically detected.\n\n"

seqbuddy_blurb += "The %s available tools:" % len(br.sb_flags)
sb_flags = OrderedDict(sorted(br.sb_flags.items(), key=lambda x: x[0]))
for func in sb_flags:
    func = re.sub("_", " ", func)
    seqbuddy_blurb += "\n    - %s" % func

alignbuddy_blurb = "Read, write, analyze, and manipulate alignment files in common formats including NEXUS, PHYLIP, " \
                   "and Stockholm. AlignBuddy can also generate alignments by wrapping common alignment programs " \
                   "(e.g., MAFFT, MUSCLE, and PAGAN), handling pre- and post- formatting automatically.\n\n"

alignbuddy_blurb += "The %s available tools:" % len(br.alb_flags)
alb_flags = OrderedDict(sorted(br.alb_flags.items(), key=lambda x: x[0]))
for func in alb_flags:
    func = re.sub("_", " ", func)
    alignbuddy_blurb += "\n    - %s" % func

phylobuddy_blurb = "Read, write, analyze, and manipulate phylogenetic tree files in Newick, Nexus, and XML formats.\n\n"

phylobuddy_blurb += "The %s available tools:" % len(br.pb_flags)
pb_flags = OrderedDict(sorted(br.pb_flags.items(), key=lambda x: x[0]))
for func in pb_flags:
    func = re.sub("_", " ", func)
    phylobuddy_blurb += "\n    - %s" % func

databasebuddy_blurb = "Search for and retrieve sequence records from NCBI, UniProt, and Ensembl. DatabaseBuddy is " \
                      "primarily used as a 'live shell', allowing the user to filter results before committing " \
                      "to downloading the actual sequences.\n\n"

databasebuddy_blurb += "The %s available tools:" % len(br.db_flags)
db_flags = OrderedDict(sorted(br.db_flags.items(), key=lambda x: x[0]))
for func in db_flags:
    func = re.sub("_", " ", func)
    databasebuddy_blurb += "\n    - %s" % func

if __name__ == '__main__':
    app = Installer(master=root)
    root.geometry("{0}x{1}+{2}+{3}".format(str(int(sw / 3)), str(int(sh / 2)), str(int(sw / 4)), str(int(sh / 4))))
    root.lift()
    root.call('wm', 'attributes', '.', '-topmost', True)
    root.after_idle(root.call, 'wm', 'attributes', '.', '-topmost', False)
    root.resizable(width=FALSE, height=FALSE)
    app.mainloop()
