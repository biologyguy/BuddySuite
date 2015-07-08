from tkinter import *
import collections
from functools import partial
import io

root = Tk()
sw = root.winfo_screenwidth()
sh = root.winfo_screenheight()

class Installer(Frame):
    container = []
    buddies = collections.OrderedDict()
    buddy_names = ["SeqBuddy", "PhyloBuddy", "AlignBuddy", "DBBuddy"]
    bs_logo = PhotoImage(file="BuddySuite-logo.gif")
    suite_logos = [PhotoImage(file="{0}-logo.gif".format(buddy)) for buddy in buddy_names]
    sb_funcs = ["guess_alphabet", "guess_format", "phylipi", "blast", "shuffle", "order_ids", "rna2dna", "dna2rna",
                "complement", "reverse_complement", "translate_cds", "select_frame", "translate6frames",
                "back_translate", "ave_seq_length", "concat_seqs", "clean_seq", "delete_metadata",
                "map_features_dna2prot", "map_features_prot2dna", "combine_features", "order_features_by_position",
                "order_features_alphabetically", "hash_sequence_ids", "pull_recs", "pull_random_recs",
                "pull_record_ends", "extract_range", "find_repeats", "delete_records", "delete_large",
                "delete_small", "delete_features", "delete_repeats", "rename", "purge", "bl2seq",
                "uppercase", "lowercase", "split_by_taxa", "molecular_weight", "isoelectric_point",
                "count_residues", "screw_formats", "raw_seq", "list_ids", "num_seqs", "merge", "split_file",
                "find_restriction_sites"]
    pb_funcs = []
    ab_funcs = []
    db_funcs = []
    funcs_list = [sb_funcs, pb_funcs, ab_funcs, db_funcs]


    for buddy in buddy_names:
        buddies[buddy] = False

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
        next_button = Button(padx=50, pady=20, text="Next", command=self.tool_select)
        next_button.place(x=int(sw/6)-75, y=sh/2-100)
        self.container.append(frame)
        self.container.append(next_button)

    def tool_select(self):
        self.clear_container()
        xpos = 10
        ypos = 10
        for indx, value in enumerate(self.buddies):
            func = partial(self.toggle_tool, value)
            button = Checkbutton(root, text=value, command=func)
            button.deselect()
            button.place(x=xpos, y=ypos)
            ypos += 30
            self.container.append(button)
        next_func = partial(self.next_tool, 0)
        next_button = Button(padx=50, pady=20, text="Next", command=next_func)
        next_button.place(x=int(sw/6)-75, y=sh/2-100)
        self.container.append(next_button)

    def next_tool(self, num):
        self.clear_container()

        if num>3:
            self.install_location()
        if not self.buddies[self.buddy_names[num]]:
            self.next_tool(num+1)
            return

        logo_label = Label(image=self.suite_logos[num], pady=20)
        logo_label.pack(side=TOP)
        self.container.append(logo_label)
        frame = Frame()
        self.container.append(frame)
        scrollbar = Scrollbar(master=frame)
        canvas = Canvas(master=frame, bd=0, yscrollcommand=scrollbar.set)
        func_box = Frame(canvas)
        r = 0
        for func in self.funcs_list[num]:
            cmd = partial(self.next_tool, num+1)
            button = Checkbutton(func_box, text=func+"()", command=cmd)
            button.deselect()
            button.grid(row=r, sticky=NW)
            r += 1
        func_box.pack()
        canvas.config(scrollregion=canvas.bbox(ALL))
        scrollbar.config(command=canvas.yview())
        canvas.pack(side=LEFT, expand=YES, fill=BOTH, padx=50, pady=50)
        scrollbar.pack(fill=Y)
        frame.pack(side=LEFT)
        next_func = partial(self.next_tool, num+1)
        next_button = Button(frame, text="Next", command=next_func)
        self.container.append(next_button)
        next_button.pack(side=BOTTOM)





    def install_location(self):
        raise SystemExit(0)

    def toggle_tool(self, name):
        self.buddies[name] = False if self.buddies[name] else True
        print(name)
        print(str(self.buddies[name]))

    def clear_container(self):
        for item in self.container:
            item.destroy()

root.title("BuddySuite Installer")
app = Installer(master=root)
root.geometry("{0}x{1}+{2}+{3}".format(str(int(sw/3)), str(int(sh/2)), str(int(sw/4)), str(int(sh/4))))
root.lift()
root.call('wm', 'attributes', '.', '-topmost', True)
root.after_idle(root.call, 'wm', 'attributes', '.', '-topmost', False)
root.resizable(width=FALSE, height=FALSE)
app.mainloop()
