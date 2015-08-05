#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from functools import partial
import shutil
from shutil import which, rmtree, copytree
from platform import *
from os import path
from configparser import *
import copy
from tempfile import TemporaryDirectory
import zipfile
from inspect import getsourcefile
import argparse

buddysuite_version = '1.alpha'

parser = argparse.ArgumentParser(prog="BuddySuite", description="The BuddySuite installer")

parser.add_argument('-v', '--version', action='version',
                    version='''\
BuddySuite {0} (2015)

Gnu General Public License, Version 2.0 (http://www.gnu.org/licenses/gpl.html)
This is free software; see the source for detailed copying conditions.
There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.
Questions/comments/concerns can be directed to Steve Bond, steve.bond@nih.gov'''.format(buddysuite_version))

parser.add_argument('-cmd', '--cmd_line', help="Force the installer to use command line version.",
                    action='store_true')

in_args = parser.parse_args()
home_dir = path.expanduser('~')
start_dir = os.getcwd()

# Prepare working space
temp_dir = TemporaryDirectory()

source_file = path.abspath(getsourcefile(lambda: 0)).split("/")
source_file = "/".join(source_file[:-1])
with zipfile.ZipFile(source_file) as ifile:
    ifile.extractall(temp_dir.name)

os.chdir(temp_dir.name)


if not in_args.cmd_line:
    try:
        from tkinter import *
        from tkinter import filedialog
    except ImportError:
        proceed = False
        print("Failed to build GUI. Package Tkinter was not found.")
        in_args.cmd_line = True

class BuddyInstall:
    @staticmethod
    def install_buddy_suite(user_system, options):
        print("Starting.")
        BuddyInstall.edit_profile()
        buddies_to_install = options[0]
        install_directory = options[1]
        shortcuts = options[2]
        for buddy in buddies_to_install:
            if not buddies_to_install[buddy]:
                shortcuts[buddy] = []

        paths_to_delete = ["resources", "blast_binaries", "Bio"]
        files_to_delete = ["SeqBuddy.py", "AlignBuddy.py", "DatabaseBuddy.py", "PhyloBuddy.py", "MyFuncs.py",
                           "config.ini"]

        if path.exists("{0}/.buddysuite/__pycache__".format(home_dir)):
            rmtree("{0}/.buddysuite/__pycache__".format(home_dir))

        all_false = True
        for buddy in buddies_to_install:
            if buddies_to_install[buddy]:
                all_false = False
        for loc in paths_to_delete:
            if path.exists("{0}/.buddysuite/{1}".format(home_dir, loc)):
                rmtree("{0}/.buddysuite/{1}".format(home_dir, loc))
        for loc in files_to_delete:
            if path.exists("{0}/.buddysuite/{1}".format(home_dir, loc)):
                os.remove("{0}/.buddysuite/{1}".format(home_dir, loc))

        if path.exists("{0}/.buddysuite".format(home_dir)) and len(os.listdir("{0}/.buddysuite".format(home_dir))) == 0:
            rmtree(path.realpath("{0}/.buddysuite".format(home_dir)))
            os.remove("{0}/.buddysuite".format(home_dir))

        myfuncs_path = "./MyFuncs.py"
        biopython_path = "./Bio"
        if not all_false:
            print("Install path: " + install_directory)
            os.makedirs(install_directory, exist_ok=True)
            print("Directory added: {0}".format(install_directory))
            if not path.exists("{0}/.buddysuite".format(home_dir)):
                os.symlink(install_directory, "{0}/.buddysuite".format(home_dir))
                print("Shortcut added: {0} ==> {1}/.buddysuite".format(install_directory, home_dir))
            if user_system in ['Darwin', 'Linux', 'Unix']:
                user_system = 'Linux' if user_system == "Unix" else user_system
                user_system = 'Win32' if user_system == "Windows" else user_system

                shutil.copy(myfuncs_path, "{0}/MyFuncs.py".format(install_directory))
                print("File added: {0}/MyFuncs.py".format(install_directory))
                copytree(biopython_path, "{0}/Bio".format(install_directory))
                print("Directory added: {0}/Bio".format(install_directory))

                binaries = ["blastn", "blastp", "blastdbcmd"] if user_system != "Win32" else \
                    ["blastn.exe", "blastp.exe", "blastdbcmd.exe"]
                for binary in binaries:
                    if not which(binary):
                        shutil.copy("blast_binaries/{0}_{1}".format(user_system, binary),
                                    "{0}/{1}".format(install_directory, binary))
                        os.chmod("{0}/{1}".format(install_directory, binary), 0o755)
                        print("{0} ==> Installed".format(binary))
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
        writer['DEFAULT'] = {'selected': {'SeqBuddy': True, 'AlignBuddy': True, 'PhyloBuddy': True, 'DatabaseBuddy': True},
                             'Install_path': {'path': '{0}/BuddySuite'.format(home_dir)},
                             'shortcuts': {'SeqBuddy': 'sb\tseqbuddy', 'AlignBuddy': 'alb\talignbuddy',
                                           'PhyloBuddy': 'pb\tphylobuddy', 'DatabaseBuddy': 'db\tDatabaseBuddy'}}

        for buddy in options[0]:
            if options[0][buddy]:
                writer['selected'][buddy] = 'True'
            else:
                writer['selected'][buddy] = 'False'

        writer["Install_path"]['path'] = options[1]

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
                       reader.get('Install_path', 'path'), {}]

            for buddy in options[0]:
                if reader['selected'][buddy] == 'True':
                    options[0][buddy] = True
                if reader['shortcuts'][buddy] != "None":
                    sc = reader['shortcuts'][buddy].split("\t")
                    options[2][buddy] = sc
                else:
                    options[2][buddy] = []

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
            with open("{0}/.profile".format(home_dir)) as file:
                contents = file.read()
                if re.search(regex, contents) is None:
                    file.close()
                    with open("{0}/.bashrc".format(home_dir), 'a') as file_write:
                        file_write.write("\n# added by BuddySuite installer\n")
                        file_write.write('export PATH=$PATH:{0}/.buddysuite\n'.format(home_dir))


def cmd_install():  # ToDo: Beef this up some, allowing as much flexibility as the graphical installer
    def ask(_prompt):
        _response = input(_prompt)
        while True:
            if _response.lower() in ["yes", "y", '']:
                return True
            elif _response.lower() in ["no", "n", "abort"]:
                return False
            else:
                print("Response not understood. Try again.")
                _response = input(_prompt)

    config = BuddyInstall.read_config_file()
    old_install_dir = None
    already_installed = None
    old_shortcuts = None
    if config is not None:
        already_installed = copy.deepcopy(config[0])
        old_install_dir = copy.deepcopy(config[1])
        install_dir = old_install_dir
        old_shortcuts = copy.deepcopy(config[2])
    buddies_to_install = {"SeqBuddy": False, "PhyloBuddy": False, "AlignBuddy": False, "DatabaseBuddy": False}
    shortcuts = {"SeqBuddy": [], "PhyloBuddy": [], "AlignBuddy": [], "DatabaseBuddy": []}

    def print_shortcuts():
        print("\tSelected Shortcuts: ")
        for buddy in shortcuts:
            if buddies_to_install[buddy]:
                print("\t{0}:\t{1}".format(buddy, str(shortcuts[buddy]).strip('[]')))

    def add_shortcut(_buddy, _shortcut):
        if buddies_to_install[_buddy] is False:
            print("Warning: {0} is not being installed. Shortcut '{1}' could not be added.".format(_buddy, _shortcut))
            return
        if which(_shortcut) is None:
            for buddy in shortcuts:
                for shortcut in shortcuts[buddy]:
                    if _shortcut == shortcut:
                        print("Warning: '{0}' is already a shortcut!".format(_shortcut))
                        return
            shortcuts[_buddy].append(_shortcut)
            print("Shortcut added: {0}\t==>\t'{1}'".format(_buddy, _shortcut))
        else:
            print("Warning: '{0}' is already in use by another program.".format(_shortcut))

    def remove_shortcut(_shortcut):
        for _buddy in shortcuts:
            if _shortcut in shortcuts[_buddy]:
                shortcuts[_buddy].remove(_shortcut)
                return True
        return False

    print("Basic terminal installer called.")
    if not ask("Would you like to proceed? ('[yes]/no') "):
        print("Installation aborted.")
        exit()

    print("Before continuing, please review our license at: \nhttp://www.gnu.org/licenses/gpl-3.0.en.html")
    if not ask("Do you accept these terms? ('[yes]/no') "):
        print("Installation aborted.")
        exit()

    if config is not None:
        print("We have detected a previous installation. Some settings will be imported.")

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
            if ask("Would you like to keep your old shortcuts? ('[yes]/no') "):
                shortcuts = old_shortcuts
            print("\tImported Shortcuts: ")
            for buddy in shortcuts:
                if buddies_to_install[buddy]:
                    print("\t{0}:\t{1}".format(buddy, str(shortcuts[buddy]).strip('[]')))

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

            if ask("Would you like to add the default shortcuts? ('[yes]/no') "):
                for buddy in default_shortcuts:
                    if buddies_to_install[buddy]:
                        for shortcut in default_shortcuts[buddy]:
                            add_shortcut(buddy, shortcut)

        print_shortcuts()

        if ask("Would you like to add or remove shortcuts? ('[yes]/no') "):
            prompt = True
            while prompt:
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
                prompt = ask("Would you like to add/remove another shortcut? ('[yes]/no') ")

        if old_install_dir is None:
            install_dir = input("Please specify an installation directory ('/BuddySuite' will automatically be appended) ")
            while True:
                if install_dir.lower() == 'abort':
                    print('Installation aborted.')
                    exit()
                try:
                    install_dir = re.sub('~', home_dir, install_dir)
                    if install_dir[0] != '/':
                        install_dir = "{0}/{1}".format(start_dir, install_dir)
                    os.makedirs("{0}/BuddySuite".format(install_dir), exist_ok=True)
                    install_dir = "{0}/BuddySuite".format(install_dir)
                    break
                except PermissionError:
                    print('Insufficient privileges to write here.')
                except:
                    print('Unknown error.')
                install_dir = input("Please specify a different installation directory (or say 'abort') ")
        else:
            install_dir = old_install_dir

    print("Please verify your settings.\n")
    for buddy in buddies_to_install:
        if already_installed is not None and already_installed[buddy]:
            if buddies_to_install[buddy]:
                print("\t{0}: Modify".format(buddy))
            else:
                print("\t{0}: Uninstall".format(buddy))
        else:
            if buddies_to_install[buddy]:
                print("\t{0}: Install".format(buddy))
            else:
                print("\t{0}: Skip".format(buddy))
    if not all_false:
        print()
        print_shortcuts()
        print("\tInstallation directory: {0}\n".format(install_dir))
    if ask("Are these settings okay? ('[yes]/abort') "):
        if all_false and config is not None:
            os.remove("{0}/.buddysuite/config.ini".format(home_dir))
        if config is not None:  # Uninstall removed shortcuts
            for buddy in old_shortcuts:
                for shortcut in old_shortcuts[buddy]:
                    if os.path.exists("{0}/.buddysuite/{1}".format(home_dir, shortcut)) and shortcut not in shortcuts[buddy]:
                        os.remove("{0}/.buddysuite/{1}".format(home_dir, shortcut))
        for buddy in buddies_to_install:
            if not buddies_to_install[buddy]:
                shortcuts[buddy] = []
        BuddyInstall.install_buddy_suite(system(), [buddies_to_install, install_dir, shortcuts])
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
    sys.path.insert(0, "./")
    root.title("BuddySuite Installer")


class Installer(Frame):
    container = []
    buddies = collections.OrderedDict()
    buddy_names = ["SeqBuddy", "PhyloBuddy", "AlignBuddy", "DatabaseBuddy"]
    for buddy in buddy_names:
        buddies[buddy] = True

    uninstall = False

    bs_logo = PhotoImage(file="{0}/BuddySuite-logo.gif".format(temp_dir.name))
    id_logo = PhotoImage(file="{0}/InstallDirectory.gif".format(temp_dir.name))
    sc_logo = PhotoImage(file="{0}/ConsoleShortcuts.gif".format(temp_dir.name))
    cs_logo = PhotoImage(file="{0}/ConfirmSelection.gif".format(temp_dir.name))

    alb_logo = PhotoImage(file="{0}/AlignBuddy-logo.gif".format(temp_dir.name))
    sb_logo = PhotoImage(file="{0}/SeqBuddy-logo.gif".format(temp_dir.name))
    pb_logo = PhotoImage(file="{0}/PhyloBuddy-logo.gif".format(temp_dir.name))
    db_logo = PhotoImage(file="{0}/DBBuddy-logo.gif".format(temp_dir.name))

    suite_logos = [sb_logo, pb_logo, alb_logo, db_logo]

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

        if which("sb") is None and not buddies["SeqBuddy"]:  # if new install of given tool, re-add default shortcuts
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

    def welcome(self, debug=None):
        self.clear_container()
        title_frame = Frame()
        welcome_label = Label(title_frame, image=self.bs_logo)
        welcome_label.pack(side=TOP)
        version_label = Label(title_frame, text="Version {0}".format(buddysuite_version))
        version_label.pack(side=RIGHT)
        title_frame.pack(pady=sh / 8)
        self.container.append(title_frame)
        button_container = Frame()
        next_button = Button(button_container, width=20, pady=20, text="Install", command=self.license)
        uninstall_button = Button(button_container, width=20, pady=20, text="Uninstall",
                                  command=self.uninstall_all)
        if self.modifying:
            next_button.config(text="Modify Installation")
            uninstall_button.pack(side=BOTTOM)
        next_button.pack(side=TOP)
        button_container.pack(side=BOTTOM, pady=40)
        self.container.append(button_container)

    def license(self):
        self.uninstall = False
        self.clear_container()
        frame = Frame(pady=75)
        scrollbar = Scrollbar(master=frame)
        license_file = open("{0}/LICENSE".format(temp_dir.name))
        license_box = Text(master=frame, wrap=WORD, yscrollcommand=scrollbar.set)
        license_box.insert(END, license_file.read())
        license_box.config(state=DISABLED)
        scrollbar.config(command=license_box.yview())
        frame.pack(side=TOP)
        license_box.pack(side=LEFT)
        scrollbar.pack(side=RIGHT, fill=Y)
        button_frame = Frame()
        next_button = Button(button_frame, padx=50, pady=20, text="I agree", command=self.next_tool)
        next_button.pack(side=TOP)
        back_button = Label(button_frame, padx=50, pady=20, text="Cancel", fg='blue', font=('Helvetica', 12,
                                                                                            'underline'))
        back_button.bind("<Button-1>", self.welcome)
        back_button.pack(side=BOTTOM)
        button_frame.pack(side=BOTTOM, pady=20)
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

        logo_label = Label(image=self.suite_logos[num], pady=20)
        logo_label.pack(side=TOP)
        self.container.append(logo_label)
        mega_frame = Frame(pady=20)
        self.container.append(mega_frame)
        frame = Frame(mega_frame, padx=50, pady=10)
        scrollbar = Scrollbar(master=frame)
        description_file = open("{0}/LICENSE".format(temp_dir.name))
        description_box = Text(master=frame, wrap=WORD, yscrollcommand=scrollbar.set)
        description_box.insert(END, description_file.read())
        description_box.config(state=DISABLED)
        scrollbar.config(command=description_box.yview())
        description_box.pack(side=LEFT)
        scrollbar.pack(side=RIGHT, fill=Y)
        button_frame = Frame(mega_frame)
        next_func = partial(self.next_tool, num + 1)
        next_button = Button(button_frame, padx=50, pady=20, text="Next", command=next_func)
        self.container.append(next_button)
        next_button.pack(side=RIGHT)
        back_func = partial(self.next_tool, num - 1)
        back_button = Button(button_frame, padx=50, pady=20, text="Back", command=back_func)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM)
        func = partial(self.toggle_tool, self.buddy_names[num])
        if not self.modifying or self.config[0][self.buddy_names[num]] is False:
            tool_button = Checkbutton(mega_frame, text="Install {0}".format(self.buddy_names[num]), command=func,
                                      pady=20)
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
                                 pady=20)
            uninstall = Radiobutton(radio_frame, text="Uninstall", value=0, variable=var, command=radiobutton_func,
                                    pady=20)
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
        label1 = Label(labelframe, text="You're currently not installing anything!", pady=5, font=('Helvetica', 16))
        label2 = Label(labelframe, text="Please go back and select at least one tool.", pady=5, font=('Helvetica', 16))
        label1.pack(side=TOP)
        label2.pack(side=BOTTOM)
        labelframe.pack(side=TOP, pady=120)
        back_func = partial(self.next_tool, 3)
        back_button = Button(root, padx=50, pady=20, text="Back", command=back_func)
        back_button.pack(side=BOTTOM, pady=40)
        self.container.append(labelframe)
        self.container.append(back_button)

    def install_location(self):
        self.clear_container()
        logo_label = Label(image=self.id_logo, pady=20)
        logo_label.pack(side=TOP)
        self.container.append(logo_label)
        frame = Frame(root, padx=50, pady=50)
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
        frame.pack(side=TOP, padx=10, pady=10, fill=BOTH)

        toggle_func = partial(self.default_directory, directory_text, browse_button)
        toggle_default = Checkbutton(frame, text="Default directory", pady=10, command=toggle_func)
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
            warning.pack(side=TOP, pady=50)
        button_frame = Frame(lower_box)
        next_func = partial(self.install_shortcuts, directory_text)
        next_button = Button(button_frame, padx=50, pady=20, text="Next", command=next_func)
        next_button.pack(side=RIGHT)
        back_func = partial(self.next_tool, 3, directory_text)
        back_button = Button(button_frame, padx=50, pady=20, text="Back", command=back_func)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM, pady=20)
        lower_box.pack(side=BOTTOM)
        self.container.append(lower_box)

    def install_shortcuts(self, in_dir=None):
        if in_dir is not None:
            self.install_dir = in_dir.get()
        self.clear_container()
        logo_label = Label(image=self.sc_logo)
        logo_label.pack(side=TOP)
        self.container.append(logo_label)

        button_frame = Frame()
        next_button = Button(button_frame, padx=50, pady=20, text="Next", command=self.confirmation)
        next_button.pack(side=RIGHT)
        back_button = Button(button_frame, padx=50, pady=20, text="Back", command=self.install_location)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM, pady=40)
        self.container.append(button_frame)

        frame = Frame()
        self.container.append(frame)
        scrollbox_frame = Frame(frame)
        scrollbar = Scrollbar(master=scrollbox_frame)
        shortcut_box = Listbox(master=scrollbox_frame, yscrollcommand=scrollbar.set, bd=2, relief=SUNKEN)
        for buddy in self.shortcuts:
            if self.buddies[buddy]:
                for shortcut in self.shortcuts[buddy]:
                    shortcut_box.insert(END, "{0} ==> {1}".format(buddy, shortcut))
        scrollbar.config(command=shortcut_box.yview())
        shortcut_box.pack(side=LEFT, fill=BOTH, expand=1)
        scrollbar.pack(side=RIGHT, fill=Y)
        scrollbox_frame.pack(side=TOP, fill=BOTH, expand=1)
        space = Frame()
        space.pack(pady=25)
        frame.pack(padx=100, expand=1, fill=BOTH, side=BOTTOM)
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
        shortcut_entry = Entry(entry_frame)
        add_func = partial(self.add_shortcut, curr_buddy, shortcut_box, shortcut_entry, debug)
        add_button = Button(entry_button_frame, text="Add", command=add_func)
        add_button.pack(side=LEFT)
        rmv_func = partial(self.remove_shortcut, curr_buddy, shortcut_box, shortcut_entry, debug)
        rmv_button = Button(entry_button_frame, text="Remove", command=rmv_func)
        rmv_button.pack(side=RIGHT)
        shortcut_entry.pack(side=LEFT, fill=X, expand=1)
        entry_button_frame.pack(side=LEFT)
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

    def add_shortcut(self, buddy, listbox, entry, debug, event=None):
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

    def remove_shortcut(self, buddy, listbox, entry, debug, event=None):
        debug.config(text='')
        text = re.sub("[^a-zA-Z0-9]", '', entry.get())
        entry.delete(0, END)
        entry.insert(END, text)
        lst = listbox.get(0, END)
        if "{0} ==> {1}".format(buddy.get(), text) in lst:
            index = lst.index("{0} ==> {1}".format(buddy.get(), text))
            listbox.delete(index)
            self.shortcuts[buddy.get()].remove(text)

    def confirmation(self):
        self.clear_container()
        logo_label = Label(image=self.cs_logo, pady=20)
        logo_label.pack(side=TOP)
        self.container.append(logo_label)
        info_frame = LabelFrame(text="Selections", bd=2, relief=SUNKEN, padx=10, pady=10)
        self.container.append(info_frame)
        os_label = Label(info_frame, text="{:<28}{}".format("Operating System:", self.user_os), font=('Courier', 13))
        dir_label = Label(info_frame, text="{:<28}{}".format("Install Directory:", self.install_dir), font=('Courier', 13))

        if not self.modifying or self.config[0]["SeqBuddy"] is False:
            install_state = "Install" if self.buddies["SeqBuddy"] else "Skip"
        else:
            install_state = "Modify" if self.buddies["SeqBuddy"] else "Uninstall"
        sb_label = Label(info_frame, text="{:<28}{}".format("SeqBuddy:", install_state), font=('Courier', 13))

        if not self.modifying or self.config[0]["PhyloBuddy"] is False:
            install_state = "Install" if self.buddies["PhyloBuddy"] else "Skip"
        else:
            install_state = "Modify" if self.buddies["PhyloBuddy"] else "Uninstall"
        pb_label = Label(info_frame, text="{:<28}{}".format("PhyloBuddy:", install_state), font=('Courier', 13))

        if not self.modifying or self.config[0]["AlignBuddy"] is False:
            install_state = "Install" if self.buddies["AlignBuddy"] else "Skip"
        else:
            install_state = "Modify" if self.buddies["AlignBuddy"] else "Uninstall"
        ab_label = Label(info_frame, text="{:<28}{}".format("AlignBuddy:", install_state), font=('Courier', 13))

        if not self.modifying or self.config[0]["DatabaseBuddy"] is False:
            install_state = "Install" if self.buddies["DatabaseBuddy"] else "Skip"
        else:
            install_state = "Modify" if self.buddies["DatabaseBuddy"] else "Uninstall"
        db_label = Label(info_frame, text="{:<28}{}".format("DatabaseBuddy:", install_state), font=('Courier', 13))

        cs_sb_label = Label(info_frame, font=('Courier', 13),
                            text="{:<28}{}".format("SeqBuddy Shortcut(s):", ", ".join(self.shortcuts["SeqBuddy"])))
        cs_ab_label = Label(info_frame, font=('Courier', 13),
                            text="{:<28}{}".format("AlignBuddy Shortcut(s):", ", ".join(self.shortcuts["AlignBuddy"])))
        cs_pb_label = Label(info_frame, font=('Courier', 13),
                            text="{:<28}{}".format("PhyloBuddy Shortcut(s):", ", ".join(self.shortcuts["PhyloBuddy"])))
        cs_db_label = Label(info_frame, font=('Courier', 13),
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
        info_frame.pack(side=TOP, anchor=NW, padx=50, pady=50, fill=BOTH)

        all_false = True
        for buddy in self.buddies:
            if self.buddies[buddy]:
                all_false = False
        if self.config is not None and all_false:
            back_func = partial(self.next_tool, 3)
        else:
            back_func = self.install_shortcuts

        button_frame = Frame()
        next_button = Button(button_frame, padx=50, pady=20, text="Install", command=self.install)
        if all_false:
            next_button.config(text="Uninstall")
        next_button.pack(side=RIGHT)
        back_button = Button(button_frame, padx=50, pady=20, text="Back", command=back_func)
        if self.uninstall:
            back_button.config(command=self.welcome)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM, pady=40)
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
                    if os.path.exists("{0}/.buddysuite/{1}".format(home_dir, shortcut)) and shortcut not in self.shortcuts[buddy]:
                        os.remove("{0}/.buddysuite/{1}".format(home_dir, shortcut))
        for buddy in self.buddies:
            if not self.buddies[buddy]:
                self.shortcuts[buddy] = []
        BuddyInstall.install_buddy_suite(self.user_system, [self.buddies, self.install_dir, self.shortcuts])
        exit()

    def clear_container(self):
        for item in self.container:
            item.destroy()

    @staticmethod
    def choose_directory(textbox):
        name = filedialog.askdirectory(parent=root, title="Select Directory", initialdir=textbox.get())
        if name is not "":
            textbox.delete(0, len(textbox.get()))
            textbox.insert(END, name+"/BuddySuite")

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


if __name__ == '__main__':
    app = Installer(master=root)
    root.geometry("{0}x{1}+{2}+{3}".format(str(int(sw/3)), str(int(sh/2)), str(int(sw/4)), str(int(sh/4))))
    root.lift()
    root.call('wm', 'attributes', '.', '-topmost', True)
    root.after_idle(root.call, 'wm', 'attributes', '.', '-topmost', False)
    root.resizable(width=FALSE, height=FALSE)
    app.mainloop()
