#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from tkinter import *
from tkinter import filedialog
import collections
from functools import partial
from shutil import which
from platform import *
from os import *
from configparser import *

root = Tk()
sw = root.winfo_screenwidth()
sh = root.winfo_screenheight()

class Installer(Frame):
    container = []
    buddies = collections.OrderedDict()
    buddy_names = ["SeqBuddy", "PhyloBuddy", "AlignBuddy", "DBBuddy"]
    for buddy in buddy_names:
        buddies[buddy] = True

    bs_logo = PhotoImage(file="./resources/images/BuddySuite-logo.gif")
    id_logo = PhotoImage(file="./resources/images/InstallDirectory.gif")
    cs_logo = PhotoImage(file="./resources/images/ConfirmSelection.gif")
    suite_logos = [PhotoImage(file="./resources/images/{0}-logo.gif".format(buddy)) for buddy in buddy_names]

    install_dir = "/usr/local/bin/BuddySuite"
    default_dir = "/usr/local/bin/BuddySuite"
    default = True
    shortcuts = {"SeqBuddy": [], "PhyloBuddy": [], "AlignBuddy": [], "DBBuddy": []}
    user_system = system()
    user_os = platform()
    print("Operating System: {0}".format(user_os))

    if BuddyInstall.read_config_file() is not None:
        config = BuddyInstall.read_config_file()
        buddies = config[0]
        install_dir = config[1]
        shortcuts = config[2]

    conflict = False
    if (which('sb') is not None) or (which('pb') is not None) or (which('alb') is not None) or (which('db') is not None):
        conflict = True
    install_shortcuts = False if conflict else True

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
        frame = Frame(padx=50, pady=100)
        scrollbar = Scrollbar(master=frame)
        license_file = open("LICENSE")
        license_box = Text(master=frame, wrap=WORD, yscrollcommand=scrollbar.set)
        license_box.insert(END, license_file.read())
        license_box.config(state=DISABLED)
        scrollbar.config(command=license_box.yview())
        frame.pack(side=TOP)
        license_box.pack(side=LEFT)
        scrollbar.pack(side=RIGHT, fill=Y)
        next_button = Button(padx=50, pady=20, text="I agree", command=self.next_tool)
        next_button.place(x=int(sw/6)-75, y=sh/2-100)
        self.container.append(frame)
        self.container.append(next_button)

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
        toggle_shortcuts = Checkbutton(text="Install Console Shortcuts", pady=10,
                                       command=self.toggle_console_shortcuts)
        if self.default:
            browse_button.config(state=DISABLED)
            directory_text.config(state=DISABLED)
            toggle_default.select()
        else:
            toggle_default.deselect()
        if self.install_shortcuts:
            toggle_shortcuts.select()
        else:
            toggle_shortcuts.deselect()
        if self.conflict:
            toggle_shortcuts.config(state=DISABLED, text="Install Console Shortcuts (ERROR: Naming conflict)")

        toggle_default.pack(side=LEFT)
        toggle_shortcuts.pack(padx=60, anchor=NW)
        self.container.append(toggle_shortcuts)
        self.container.append(frame)
        button_frame = Frame()
        next_func = partial(self.confirmation, directory_text)
        next_button = Button(button_frame, padx=50, pady=20, text="Next", command=next_func)
        next_button.pack(side=RIGHT)
        back_func = partial(self.next_tool, 3, directory_text)
        back_button = Button(button_frame, padx=50, pady=20, text="Back", command=back_func)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM, pady=40)
        self.container.append(button_frame)

    def toggle_tool(self, name):
        self.buddies[name] = False if self.buddies[name] else True
        print(name)
        print(str(self.buddies[name]))

    def toggle_console_shortcuts(self):
        self.install_shortcuts = False if self.install_shortcuts else True

    def confirmation(self, in_dir):
        self.install_dir = in_dir.get()
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
        db_label = Label(info_frame, text="Install DBBuddy: {0}".format(self.buddies["DBBuddy"]))
        dir_label = Label(info_frame, text="Install Directory: {0}".format(self.install_dir))
        if self.conflict:
            short = "Naming conflict"
        else:
            short = self.install_shortcuts
        cs_label = Label(info_frame, text="Install Console Shortcuts: {0}".format(short))
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
        back_button = Button(button_frame, padx=50, pady=20, text="Back", command=self.install_location)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM, pady=40)
        self.container.append(button_frame)

    def install(self):
        BuddyInstall.make_config_file([self.buddies, self.install_dir, self.shortcuts])
        raise SystemExit(0)

    def clear_container(self):
        for item in self.container:
            item.destroy()

    def choose_directory(self, textbox):
        name = filedialog.askdirectory(parent=root, title="Select Directory", initialdir=textbox.get())
        if name is not "":
            textbox.delete(0, len(textbox.get()))
            textbox.insert(END, name)

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


class BuddyInstall:
    sys.path.insert(0, "./")

    @staticmethod
    def install_buddy_suite(user_system, options):
        buddies_to_install = options[0]
        install_directory = options[1]
        shortcuts = options[2]

        myfuncs_path = "./MyFuncs.py"
        biopython_path = "./Bio"
        blast_path = "./blast binaries"
        resource_path = "./resources"
        if user_system in ['Darwin', 'Linux', 'Unix']:
            mkdir(install_directory)
            copytree(biopython_path, "{0}/".format(install_directory))
            copytree(myfuncs_path, "{0}/".format(install_directory))
            copytree(blast_path, "{0}/".format(install_directory))
            copytree(resource_path, "{0}/".format(install_directory))
            for buddy in buddies_to_install:
                if buddies_to_install[buddy]:
                    copytree("./{0}.py".format(buddy), "{0}/{1}.py".format(install_directory, buddy))
                    for shortcut in shortcuts[buddy]:
                        if which(shortcut) is None:
                            os.symlink("{0}/{1}.py".format(install_directory, buddy),
                                   "/usr/local/bin/{0}".format(shortcut))
        elif user_system == 'Windows':
            return
        BuddyInstall.make_config_file(options)

    @staticmethod
    def make_config_file(options):
        writer = ConfigParser()
        writer['DEFAULT'] = {'selected': {'SeqBuddy': True, 'AlignBuddy': True, 'PhyloBuddy': True, 'DBBuddy': True},
                             'Install_path': {'path': '/usr/local/bin/BuddySuite'},
                             'shortcuts': {'SeqBuddy': 'sb\n\tseqbuddy', 'AlignBuddy': 'alb\n\talignbuddy',
                                           'PhyloBuddy': 'pb\n\tphylobuddy', 'DBBuddy': 'db\n\tdbbuddy'}}

        for buddy in options[0]:
            if options[0][buddy]:
                writer['selected'][buddy] = True
            else:
                writer['selected'][buddy] = False

        writer["Install_path"]['path'] = options[1]

        for buddy in options[2]:
            sc = ''
            for shortcut in options[2][buddy]:
                sc += shortcut + "\n\t"
            writer["{0}_shortcuts".format(buddy)] = sc if sc != '' else 'None'

        with open("{0}/resources/config.ini".format(options[0]), 'w') as configfile:
            writer.write(configfile)

    @staticmethod
    def read_config_file():
        if path.exists("./resources/config.ini"):
            reader = ConfigParser()
            reader.read_string("./resources/config.ini")
            options = [{"SeqBuddy": False, "AlignBuddy": False, "PhyloBuddy": False, "DBBuddy": False},
                       reader['Install_path']['path'], {}]

            for buddy in options[0]:
                if reader['selected'][buddy] == 'True':
                    options[0][buddy] = True
                if reader['shortcuts'][buddy] != "None":
                    sc = reader['shortcuts'][buddy].split("\n\t")
                    options[2][buddy] = sc
                else:
                    options[2][buddy] = []

            return options

        else:
            return None

root.title("BuddySuite Installer")
app = Installer(master=root)
root.geometry("{0}x{1}+{2}+{3}".format(str(int(sw/3)), str(int(sh/2)), str(int(sw/4)), str(int(sh/4))))
root.lift()
root.call('wm', 'attributes', '.', '-topmost', True)
root.after_idle(root.call, 'wm', 'attributes', '.', '-topmost', False)
root.resizable(width=FALSE, height=FALSE)
app.mainloop()

