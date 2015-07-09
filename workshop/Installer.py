from tkinter import *
from tkinter import filedialog
import collections
from functools import partial

root = Tk()
sw = root.winfo_screenwidth()
sh = root.winfo_screenheight()

class Installer(Frame):
    container = []
    buddies = collections.OrderedDict()
    buddy_names = ["SeqBuddy", "PhyloBuddy", "AlignBuddy", "DBBuddy"]
    for buddy in buddy_names:
        buddies[buddy] = True

    bs_logo = PhotoImage(file="BuddySuite-logo.gif")
    id_logo = PhotoImage(file="InstallDirectory.gif")
    cs_logo = PhotoImage(file="ConfirmSelection.gif")
    suite_logos = [PhotoImage(file="{0}-logo.gif".format(buddy)) for buddy in buddy_names]

    install_dir = "/usr/local/bin/BuddySuite"
    default_dir = "/usr/local/bin/BuddySuite"
    default = True

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
        if self.default:
            browse_button.config(state=DISABLED)
            directory_text.config(state=DISABLED)
            toggle_default.select()
        else:
            toggle_default.deselect()
        toggle_default.pack(side=LEFT)
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

    def confirmation(self, in_dir):
        self.install_dir = in_dir.get()
        self.clear_container()
        logo_label = Label(image=self.cs_logo, pady=20)
        logo_label.pack(side=TOP)
        self.container.append(logo_label)
        info_frame = Frame()
        self.container.append(info_frame)
        sb_label = Label(info_frame, text="Install SeqBuddy: {0}".format(self.buddies["SeqBuddy"]))
        pb_label = Label(info_frame, text="Install PhyloBuddy: {0}".format(self.buddies["PhyloBuddy"]))
        ab_label = Label(info_frame, text="Install AlignBuddy: {0}".format(self.buddies["AlignBuddy"]))
        db_label = Label(info_frame, text="Install DBBuddy: {0}".format(self.buddies["DBBuddy"]))
        dir_label = Label(info_frame, text="Install Directory: {0}".format(self.install_dir))
        sb_label.grid(row=0, sticky=NW)
        pb_label.grid(row=1, sticky=NW)
        ab_label.grid(row=2, sticky=NW)
        db_label.grid(row=3, sticky=NW)
        dir_label.grid(row=4, sticky=NW)
        info_frame.pack(side=TOP, anchor=NW, padx=50, pady=50)
        button_frame = Frame()
        next_button = Button(button_frame, padx=50, pady=20, text="Install", command=self.install)
        next_button.pack(side=RIGHT)
        back_button = Button(button_frame, padx=50, pady=20, text="Back", command=self.install_location)
        back_button.pack(side=LEFT)
        button_frame.pack(side=BOTTOM, pady=40)
        self.container.append(button_frame)

    def install(self):
        # run install script
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


root.title("BuddySuite Installer")
app = Installer(master=root)
root.geometry("{0}x{1}+{2}+{3}".format(str(int(sw/3)), str(int(sh/2)), str(int(sw/4)), str(int(sh/4))))
root.lift()
root.call('wm', 'attributes', '.', '-topmost', True)
root.after_idle(root.call, 'wm', 'attributes', '.', '-topmost', False)
root.resizable(width=FALSE, height=FALSE)
app.mainloop()

