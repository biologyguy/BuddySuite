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

try:
    os.makedirs("/usr/bin/temp/asfkdsgeriugengdfsvjkdvjlirutghjdfnb")
    shutil.rmtree("/usr/bin/temp/asfkdsgeriugengdfsvjkdvjlirutghjdfnb")
except PermissionError:
    print("Error: You need to run the program as a superuser/administrator.")
    raise SystemExit

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
            raise SystemExit
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
            raise SystemExit
        else:
            response = input("Response not understood. Try again. \nWould you like to proceed? ('yes/no')")

class BuddyInstall:

    @staticmethod
    def install_buddy_suite(user_system, options=None):
        if options is not None:
            buddies_to_install = options[0]
            install_directory = options[1]
            shortcuts = options[2]
        else:
            buddies_to_install = {"SeqBuddy": True, "AlignBuddy": True, "PhyloBuddy": True, "DatabaseBuddy": True}
            install_directory = "/usr/local/bin/BuddySuite"
            shortcuts = {"SeqBuddy": [], "AlignBuddy": [], "PhyloBuddy": [], "DatabaseBuddy": []}
            if which("sb") is not None:
                shortcuts["SeqBuddy"] = ["sb"]
            if which("alb") is not None:
                shortcuts["AlignBuddy"] = ["alb"]
            if which("pb") is not None:
                shortcuts["PhyloBuddy"] = ["pb"]
            if which("db") is not None:
                shortcuts["DatabaseBuddy"] = ["db"]

        paths_to_delete = ["/resources", "blast_binaries", "Bio"]
        files_to_delete = ["SeqBuddy.py", "AlignBuddy.py", "DatabaseBuddy.py", "PhyloBuddy.py", "MyFuncs.py"]
        for loc in paths_to_delete:
            if path.exists("/usr/local/bin/buddysuite/{0}".format(loc)):
                shutil.rmtree("/usr/local/bin/buddysuite/{0}".format(loc))
        for loc in files_to_delete:
            if path.exists("/usr/local/bin/buddysuite/{0}".format(loc)):
                os.remove("/usr/local/bin/buddysuite/{0}".format(loc))


        myfuncs_path = "./MyFuncs.py"
        biopython_path = "./Bio"
        blast_path = "./blast_binaries"
        resource_path = "./resources"
        print(install_directory)
        if not path.exists(install_directory):
            mkdir(install_directory)
        if user_system in ['Darwin', 'Linux', 'Unix']:
            shutil.copy(myfuncs_path, "{0}/MyFuncs.py".format(install_directory))
            BuddyInstall.copytree(resource_path, "{0}/resources".format(install_directory))
            BuddyInstall.copytree(biopython_path, "{0}/Bio".format(install_directory))
            BuddyInstall.copytree(blast_path, "{0}/blast_binaries".format(install_directory))
            for buddy in buddies_to_install:
                if buddies_to_install[buddy]:
                    shutil.copy("./{0}.py".format(buddy), "{0}/{1}.py".format(install_directory, buddy))
                    for shortcut in shortcuts[buddy]:
                        if which(shortcut) is None:
                            os.symlink("{0}/{1}.py".format(install_directory, buddy),
                                       "/usr/local/bin/{0}".format(shortcut))
            if not path.exists("/usr/local/bin/buddysuite/"):
                os.symlink(install_directory, "/usr/local/bin/buddysuite")

        elif user_system == 'Windows':
            return

        BuddyInstall.make_config_file(options)

    @staticmethod
    def make_config_file(options):
        writer = ConfigParser()
        writer.add_section('selected')
        writer.add_section('Install_path')
        writer.add_section('shortcuts')
        writer['DEFAULT'] = {'selected': {'SeqBuddy': True, 'AlignBuddy': True, 'PhyloBuddy': True, 'DatabaseBuddy': True},
                             'Install_path': {'path': '/usr/local/bin/.BuddySuite'},
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

    @staticmethod
    def read_config_file():
        if path.exists("/usr/local/bin/buddysuite/config.ini"):
            reader = ConfigParser()
            reader.read("/usr/local/bin/buddysuite/config.ini")

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

if not hard_install:
    root = Tk()
    sw = root.winfo_screenwidth()
    sh = root.winfo_screenheight()
    sys.path.insert(0, "./")
    root.title("BuddySuite Installer")
else:
    BuddyInstall.install_buddy_suite(system())
    raise SystemExit

class Installer(Frame):
    container = []
    buddies = collections.OrderedDict()
    buddy_names = ["SeqBuddy", "PhyloBuddy", "AlignBuddy", "DatabaseBuddy"]
    for buddy in buddy_names:
        buddies[buddy] = True

    bs_logo = PhotoImage(file="./resources/images/BuddySuite-logo.gif")
    id_logo = PhotoImage(file="./resources/images/InstallDirectory.gif")
    sc_logo = PhotoImage(file="./resources/images/ConsoleShortcuts.gif")
    cs_logo = PhotoImage(file="./resources/images/ConfirmSelection.gif")
    suite_logos = [PhotoImage(file="./resources/images/{0}-logo.gif".format(buddy)) for buddy in buddy_names]

    install_dir = "/usr/local/bin/.BuddySuite"
    default_dir = "/usr/local/bin/.BuddySuite"
    default = True
    shortcuts = {"SeqBuddy": ['sb', 'seqbuddy'], "PhyloBuddy": ['pb', 'phylobuddy'],
                 "AlignBuddy": ['alb', 'alignbuddy'], "DatabaseBuddy": ['db', 'dbbuddy']}
    user_system = system()
    user_os = platform()
    print("Operating System: {0}".format(user_os))

    if BuddyInstall.read_config_file() is not None:
        config = BuddyInstall.read_config_file()
        buddies = config[0]
        install_dir = config[1]
        shortcuts = config[2]
        for buddy in shortcuts:
            for shortcut in shortcuts[buddy]:
                os.remove("/usr/local/bin/{0}".format(shortcut))
    original_shortcuts = copy.deepcopy(shortcuts)
    conflict = False

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.welcome()
        self.pack()

    def welcome(self):
        welcome_label = Label(image=self.bs_logo)
        welcome_label.pack(pady=sh/8, side=TOP)
        self.container.append(welcome_label)
        next_button = Button(padx=75, pady=30, text="Install", command=self.license)
        next_button.pack(side=TOP)
        self.container.append(next_button)

    def license(self):
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
        if num<0:
            self.license()
            return
        if num>3:
            self.install_location()
            return
        logo_label = Label(image=self.suite_logos[num], pady=20)
        logo_label.pack(side=TOP)
        self.container.append(logo_label)
        mega_frame = Frame(pady=20)
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
        next_button.pack(side=RIGHT)
        back_func = partial(self.next_tool, num-1)
        back_button = Button(button_frame, padx=50, pady=20, text="Back", command=back_func)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM)
        func = partial(self.toggle_tool, self.buddy_names[num])
        tool_button = Checkbutton(mega_frame, text="Install {0}".format(self.buddy_names[num]), command=func, pady=20)
        if self.buddies[self.buddy_names[num]]:
            tool_button.select()
        else:
            tool_button.deselect()
        tool_button.pack(side=BOTTOM)
        frame.pack(side=TOP)
        mega_frame.pack(side=TOP)
        self.container.append(mega_frame)
        self.container.append(next_button)

    def toggle_tool(self, name):
        self.buddies[name] = False if self.buddies[name] else True
        print(name)
        print(str(self.buddies[name]))

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

        if self.default:
            browse_button.config(state=DISABLED)
            directory_text.config(state=DISABLED)
            toggle_default.select()
        else:
            toggle_default.deselect()

        toggle_default.pack(side=LEFT)
        self.container.append(frame)
        button_frame = Frame()
        next_func = partial(self.install_shortcuts, directory_text)
        next_button = Button(button_frame, padx=50, pady=20, text="Next", command=next_func)
        next_button.pack(side=RIGHT)
        back_func = partial(self.next_tool, 3, directory_text)
        back_button = Button(button_frame, padx=50, pady=20, text="Back", command=back_func)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM, pady=40)
        self.container.append(button_frame)

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
                    shortcut_box.insert(END,"{0} ==> {1}".format(buddy, shortcut))
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
        curr_buddy.set("SeqBuddy")
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
        cs_label = Label(info_frame, text="Console Shortcuts: {0}".format("[Placeholder]"))
        os_label.grid(row=0, sticky=NW)
        sb_label.grid(row=1, sticky=NW)
        pb_label.grid(row=2, sticky=NW)
        ab_label.grid(row=3, sticky=NW)
        db_label.grid(row=4, sticky=NW)
        dir_label.grid(row=5, sticky=NW)
        cs_label.grid(row=6, sticky=NW)
        info_frame.pack(side=TOP, anchor=NW, padx=50, pady=50, fill=BOTH)
        button_frame = Frame()
        next_button = Button(button_frame, padx=50, pady=20, text="Install", command=self.install)
        next_button.pack(side=RIGHT)
        back_button = Button(button_frame, padx=50, pady=20, text="Back", command=self.install_shortcuts)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM, pady=40)
        self.container.append(button_frame)

    def install(self):
        for buddy in self.buddies:
            if not self.buddies[buddy]:
                self.shortcuts[buddy] = []
        BuddyInstall.install_buddy_suite(self.user_system, [self.buddies, self.install_dir, self.shortcuts])
        raise SystemExit(0)

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