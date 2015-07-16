from os import *
from shutil import *
from configparser import *
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
        writer['DEFAULT'] = {'SeqBuddy': True, 'AlignBuddy': True, 'PhyloBuddy': True, 'DBBuddy': True,
                             'Install_path': '/usr/local/bin/BuddySuite', 'SeqBuddy_shortcuts': 'sb, seqbuddy',
                             'AlignBuddy_shortcuts': 'alb, alignbuddy', 'PhyloBuddy_shortcuts': 'pb, phylobuddy',
                             'DBBuddy_shortcuts': 'db, dbbuddy'}

        for buddy in options[0]:
            if options[0][buddy]:
                writer[buddy] = True
            else:
                writer[buddy] = False

        writer["Install_path"] = options[1]

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
                       reader['Install_path'], {}]

            for buddy in options[0]:
                if reader[buddy] == 'True':
                    options[0][buddy] = True
                if reader["{0}_shortcuts".format(buddy)] != "None":
                    sc = reader["{0}_shortcuts"].split("\n\t")
                    options[2][buddy] = sc
                else:
                    options[2][buddy] = []

            return options

        else:
            return None
