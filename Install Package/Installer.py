#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import collections
from functools import partial
from shutil import *
import shutil
from shutil import copy, which, copy2, copystat
from platform import *
from os import path, mkdir
from configparser import *
import copy
from re import sub
import stat

home_dir = path.expanduser('~')
hard_install = False

try:
    from tkinter import *
    from tkinter import filedialog

except ImportError:
    proceed = False
    print("Failed to build GUI. Package Tkinter was not found.")
    print("Installation will proceed through the terminal. You will not be allowed to change any options.")
    response = input("Would you like to proceed? ('yes/no')")
    while True:
        if response.lower() in ["yes", "y"]:
            break
        elif response.lower() in ["no", "n"]:
            print("Installation aborted.")
            exit()
        else:
            response = input("Response not understood. Try again. \nWould you like to proceed? ('yes/no')")
    print("Before continuing, please review our license at: \nhttp://www.gnu.org/licenses/gpl-3.0.en.html")
    response = input("Do you accept these terms? ('yes/no')")
    while True:
        if response.lower() in ["yes", "y"]:
            hard_install = True
            break
        elif response.lower() in ["no", "n"]:
            print("Installation aborted.")
            exit()
        else:
            response = input("Response not understood. Try again. \nWould you like to proceed? ('yes/no')")

class BuddyInstall:

    @staticmethod
    def install_buddy_suite(user_system, options=None):
        print("Starting installation.")
        BuddyInstall.edit_profile()
        if options is not None:
            buddies_to_install = options[0]
            install_directory = options[1]
            shortcuts = options[2]
            for buddy in buddies_to_install:
                if not buddies_to_install[buddy]:
                    shortcuts[buddy] = []

        else:
            buddies_to_install = {"SeqBuddy": True, "AlignBuddy": True, "PhyloBuddy": True, "DatabaseBuddy": True}
            install_directory = "{0}/BuddySuite".format(home_dir)

            shortcuts = {"SeqBuddy": [], "AlignBuddy": [], "PhyloBuddy": [], "DatabaseBuddy": []}
            if which("sb") is None:  # if hard install, use default shortcuts
                shortcuts["SeqBuddy"].append("sb")
            if which("seqbuddy") is None:
                shortcuts["SeqBuddy"].append("seqbuddy")
            if which("alb") is None:
                shortcuts["AlignBuddy"].append("alb")
            if which("alignbuddy") is None:
                shortcuts["AlignBuddy"].append("alignbuddy")
            if which("pb") is None:
                shortcuts["PhyloBuddy"].append("pb")
            if which("phylobuddy") is None:
                shortcuts["PhyloBuddy"].append("phylobuddy")
            if which("db") is None:
                shortcuts["DatabaseBuddy"].append("db")
            if which("databasebuddy") is None:
                shortcuts["DatabaseBuddy"].append("database")

            options = [buddies_to_install, install_directory, shortcuts]

        paths_to_delete = ["resources", "blast_binaries", "Bio"]
        files_to_delete = ["SeqBuddy.py", "AlignBuddy.py", "DatabaseBuddy.py", "PhyloBuddy.py", "MyFuncs.py",
                           "config.ini"]

        all_false = True
        for buddy in buddies_to_install:
            if buddies_to_install[buddy]:
                all_false = False
        for loc in paths_to_delete:
            if path.exists("{0}/buddysuite/{1}".format(home_dir, loc)):
                shutil.rmtree("{0}/buddysuite/{1}".format(home_dir, loc))
        for loc in files_to_delete:
            if path.exists("{0}/buddysuite/{1}".format(home_dir, loc)):
                os.remove("{0}/buddysuite/{1}".format(home_dir, loc))

        if path.exists("{0}/buddysuite".format(home_dir)) and len(os.listdir("{0}/buddysuite".format(home_dir))) == 0:
            shutil.rmtree(path.realpath("{0}/buddysuite".format(home_dir)))
            os.remove("{0}/buddysuite".format(home_dir))

        myfuncs_path = "./MyFuncs.py"
        biopython_path = "./Bio"
        blast_path = "./blast_binaries"
        resource_path = "./resources"
        if not all_false:
            print("Install path: " + install_directory)
            if not path.exists(install_directory):
                os.makedirs(install_directory)
                print("Directory added: {0}".format(install_directory))
                if not path.exists("{0}/buddysuite/"):
                    os.symlink(install_directory, "{0}/buddysuite".format(home_dir))
                    print("Shortcut added: {0} ==> {1}/buddysuite".format(install_directory, home_dir))
            if user_system in ['Darwin', 'Linux', 'Unix']:
                shutil.copy(myfuncs_path, "{0}/MyFuncs.py".format(install_directory))
                print("File added: {0}/MyFuncs.py".format(install_directory))
                BuddyInstall.copytree(resource_path, "{0}/resources".format(install_directory))
                print("Directory added: {0}/resources".format(install_directory))
                BuddyInstall.copytree(biopython_path, "{0}/Bio".format(install_directory))
                print("Directory added: {0}/Bio".format(install_directory))
                BuddyInstall.copytree(blast_path, "{0}/blast_binaries".format(install_directory))
                print("Directory added: {0}/blast_binaries".format(install_directory))
                for buddy in buddies_to_install:
                    if buddies_to_install[buddy]:
                        shutil.copy("./{0}.py".format(buddy), "{0}/{1}.py".format(install_directory, buddy))
                        print("File added: {0}/{1}.py".format(install_directory, buddy))
                        for shortcut in shortcuts[buddy]:
                            if which(shortcut) is None:
                                os.symlink("{0}/{1}.py".format(install_directory, buddy),
                                           "{0}/{1}".format(install_directory, shortcut))
                                print("{0}/{1}".format(install_directory, shortcut))
                                print("Shortcut added: {0} ==> {1}".format(buddy, shortcut))

            elif user_system == 'Windows':
                print("Windows not supported at the moment.")
                return

            BuddyInstall.make_config_file(options)
        print("Finished.")

    @staticmethod
    def uninstall_buddy_suite():
        if path.exists("{0}/buddysuite".format(home_dir)):
            paths_to_delete = ["resources", "blast_binaries", "Bio"]
            files_to_delete = ["SeqBuddy.py", "AlignBuddy.py", "DatabaseBuddy.py", "PhyloBuddy.py", "MyFuncs.py",
                               "config.ini"]
            shortcuts = BuddyInstall.read_config_file()[2]
            for buddy in shortcuts:
                for shortcut in shortcuts[buddy]:
                    if path.exists("{0}/{1}".format(home_dir, shortcut)):
                        os.remove("{0}/{1}".format(home_dir, shortcut))

            for loc in paths_to_delete:
                if path.exists("{0}/buddysuite/{1}".format(home_dir, loc)):
                    shutil.rmtree("{0}/buddysuite/{1}".format(home_dir, loc))
            for loc in files_to_delete:
                if path.exists("{0}/buddysuite/{1}".format(home_dir, loc)):
                    os.remove("{0}/buddysuite/{1}".format(home_dir, loc))

            if len(os.listdir("{0}/buddysuite".format(home_dir))) == 0:
                shutil.rmtree(path.realpath("{0}/buddysuite".format(home_dir)))
                os.remove("{0}/buddysuite".format(home_dir))

        print("BuddySuite uninstalled.")
        exit()

    @staticmethod
    def make_config_file(options):
        print("Making config file.")
        writer = ConfigParser()
        writer.add_section('selected')
        writer.add_section('Install_path')
        writer.add_section('shortcuts')
        writer['DEFAULT'] = {'selected': {'SeqBuddy': True, 'AlignBuddy': True, 'PhyloBuddy': True, 'DatabaseBuddy': True},
                             'Install_path': {'path': '{0}/.BuddySuite'.format(home_dir)},
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
        print("Config file written to "+"{0}/config.ini".format(options[1]))

    @staticmethod
    def read_config_file():
        if path.exists("{0}/buddysuite/config.ini".format(home_dir)):
            reader = ConfigParser()
            reader.read("{0}/buddysuite/config.ini".format(home_dir))

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
    def copytree(src, dst, symlinks = False, ignore=None):
        if not os.path.exists(dst):
            os.makedirs(dst)
            copystat(src, dst)
        lst = os.listdir(src)
        if ignore:
            excl = ignore(src, lst)
            lst = [x for x in lst if x not in excl]
        for item in lst:
            s = os.path.join(src, item)
            d = os.path.join(dst, item)
            if symlinks and os.path.islink(s):
                if os.path.exists(d):
                    os.remove(d)
                os.symlink(os.readlink(s), d)
                try:
                    st = os.lstat(s)
                    mode = stat.S_IMODE(st.st_mode)
                    os.lchmod(d, mode)
                except:
                    pass
            elif os.path.isdir(s):
                BuddyInstall.copytree(s, d, symlinks, ignore)
            else:
                copy2(s, d)

    @staticmethod
    def edit_profile():
        with open("{0}/.profile".format(home_dir)) as file:
            if 'export PATH="{0}/buddysuite/"'.format(home_dir) not in file.read():
                file.close()
                with open("{0}/.profile".format(home_dir), 'a') as file_write:
                    file_write.write("# added by BuddySuite installer")
                    file_write.write('export PATH="{0}/buddysuite/"'.format(home_dir))



if not hard_install:
    root = Tk()
    sw = root.winfo_screenwidth()
    sh = root.winfo_screenheight()
    sys.path.insert(0, "./")
    root.title("BuddySuite Installer")
else:
    BuddyInstall.install_buddy_suite(system())
    exit()

class Installer(Frame):
    container = []
    buddies = collections.OrderedDict()
    buddy_names = ["SeqBuddy", "PhyloBuddy", "AlignBuddy", "DatabaseBuddy"]
    for buddy in buddy_names:
        buddies[buddy] = True

    uninstall = False

    bs_logo = PhotoImage(file="./resources/images/BuddySuite-logo.gif")
    id_logo = PhotoImage(file="./resources/images/InstallDirectory.gif")
    sc_logo = PhotoImage(file="./resources/images/ConsoleShortcuts.gif")
    cs_logo = PhotoImage(file="./resources/images/ConfirmSelection.gif")
    suite_logos = [PhotoImage(file="./resources/images/{0}-logo.gif".format(buddy)) for buddy in buddy_names]

    install_dir = "{0}/.BuddySuite".format(home_dir)
    default_dir = "{0}/.BuddySuite".format(home_dir)
    default = True
    shortcuts = {"SeqBuddy": ['sb', 'seqbuddy'], "PhyloBuddy": ['pb', 'phylobuddy'],
                 "AlignBuddy": ['alb', 'alignbuddy'], "DatabaseBuddy": ['db', 'dbbuddy']}
    for buddy in shortcuts:
        for shortcut in shortcuts[buddy]:
            if which(shortcut) is not None:
                shortcuts[buddy].remove(shortcut)
    user_system = system()
    user_os = platform()

    config = None
    if BuddyInstall.read_config_file() is not None:
        config = BuddyInstall.read_config_file()
        buddies = config[0]
        install_dir = config[1]
        shortcuts = config[2]

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
            shortcuts["DatabaseBuddy"].append("database")

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

    def welcome(self):
        self.clear_container()
        welcome_label = Label(image=self.bs_logo)
        welcome_label.pack(pady=sh/8, side=TOP)
        self.container.append(welcome_label)
        button_container = Frame()
        next_button = Button(button_container, width=20, pady=20, text="Install", command=self.license)
        uninstall_button = Button(button_container, width=20, pady=20, text="Uninstall",
                                  command=self.uninstall_all)
        if self.config is not None:
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
        license_file = open("LICENSE")
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
        back_button.bind("<Button-1>", exit)
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
                if self.config is not None:
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
        description_file = open("LICENSE")
        description_box = Text(master=frame, wrap=WORD, yscrollcommand=scrollbar.set)
        description_box.insert(END, description_file.read())
        description_box.config(state=DISABLED)
        scrollbar.config(command=description_box.yview())
        description_box.pack(side=LEFT)
        scrollbar.pack(side=RIGHT, fill=Y)
        button_frame = Frame(mega_frame)
        next_func = partial(self.next_tool, num+1)
        next_button = Button(button_frame, padx=50, pady=20, text="Next", command=next_func)
        self.container.append(next_button)
        next_button.pack(side=RIGHT)
        back_func = partial(self.next_tool, num-1)
        back_button = Button(button_frame, padx=50, pady=20, text="Back", command=back_func)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM)
        func = partial(self.toggle_tool, self.buddy_names[num])
        config_file = BuddyInstall.read_config_file()
        if config_file is None or config_file[0][self.buddy_names[num]] is False:
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
        if self.config is not None:
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

    def click_shortcut(self, entry, var, event):
        event.widget.activate(event.widget.nearest(event.y))
        entry.delete(0, END)
        text = event.widget.get(ACTIVE).split(" ==> ")
        var.set(text[0])
        entry.insert(END, text[1])

    def add_shortcut(self, buddy, listbox, entry, debug):
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
        os_label = Label(info_frame, text="Operating System: {0}".format(self.user_os))
        sb_label = Label(info_frame, text="Install SeqBuddy: {0}".format(self.buddies["SeqBuddy"]))
        pb_label = Label(info_frame, text="Install PhyloBuddy: {0}".format(self.buddies["PhyloBuddy"]))
        ab_label = Label(info_frame, text="Install AlignBuddy: {0}".format(self.buddies["AlignBuddy"]))
        db_label = Label(info_frame, text="Install DatabaseBuddy: {0}".format(self.buddies["DatabaseBuddy"]))
        dir_label = Label(info_frame, text="Install Directory: {0}".format(self.install_dir))
        cs_sb_label = Label(info_frame,
                            text="SeqBuddy Shortcuts: {0}".format(str(self.shortcuts["SeqBuddy"])))
        cs_ab_label = Label(info_frame,
                            text="AlignBuddy Shortcuts: {0}".format(str(self.shortcuts["AlignBuddy"])))
        cs_pb_label = Label(info_frame,
                            text="PhyloBuddy Shortcuts: {0}".format(str(self.shortcuts["PhyloBuddy"])))
        cs_db_label = Label(info_frame,
                            text="DatabaseBuddy Shortcuts: {0}".format(str(self.shortcuts["DatabaseBuddy"])))
        os_label.grid(row=0, sticky=NW)
        sb_label.grid(row=1, sticky=NW)
        pb_label.grid(row=2, sticky=NW)
        ab_label.grid(row=3, sticky=NW)
        db_label.grid(row=4, sticky=NW)
        dir_label.grid(row=5, sticky=NW)
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
            os.remove("{0}/buddysuite/config.ini".format(home_dir))
        if self.config is not None:  # Uninstall removed shortcuts
            for buddy in self.original_shortcuts:
                for shortcut in self.original_shortcuts[buddy]:
                    if os.path.exists("{0}/buddysuite/{1}".format(home_dir, shortcut)) and shortcut not in self.shortcuts[buddy]:
                        os.remove("{0}/buddysuite/{1}".format(home_dir, shortcut))
        for buddy in self.buddies:
            if not self.buddies[buddy]:
                self.shortcuts[buddy] = []
        BuddyInstall.install_buddy_suite(self.user_system, [self.buddies, self.install_dir, self.shortcuts])
        exit()

    def clear_container(self):
        for item in self.container:
            item.destroy()

    def choose_directory(self, textbox):
        name = filedialog.askdirectory(parent=root, title="Select Directory", initialdir=textbox.get())
        if name is not "":
            textbox.delete(0, len(textbox.get()))
            textbox.insert(END, name+"/.BuddySuite")

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


app = Installer(master=root)
root.geometry("{0}x{1}+{2}+{3}".format(str(int(sw/3)), str(int(sh/2)), str(int(sw/4)), str(int(sh/4))))
root.lift()
root.call('wm', 'attributes', '.', '-topmost', True)
root.after_idle(root.call, 'wm', 'attributes', '.', '-topmost', False)
root.resizable(width=FALSE, height=FALSE)
app.mainloop()