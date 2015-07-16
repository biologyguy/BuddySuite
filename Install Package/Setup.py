from os import *
from shutil import *

sys.path.insert(0, "./")

def install_buddy_suite(user_system, options):
    buddies_to_install = options[0]
    install_directory = options[1]
    install_shortcuts = options[2]
    myfuncs_path = "./MyFuncs.py"
    biopython_path = "./Bio"
    blast_path = "./blast"
    resource_path = "./resources"
    buddy_to_shortcut = {"SeqBuddy": 'sb', "AlignBuddy": 'ab', 'PhyloBuddy': 'pb', 'DBBuddy': 'db'}
    if user_system in ['Darwin', 'Linux', 'Unix']:
        mkdir("{0}/BuddySuite".format(install_directory))
        copytree(biopython_path, "{0}/BuddySuite/".format(install_directory))
        copytree(myfuncs_path, "{0}/BuddySuite/".format(install_directory))
        copytree(blast_path, "{0}/BuddySuite/".format(install_directory))
        copytree(resource_path, "{0}/BuddySuite/".format(install_directory))
        for buddy in buddies_to_install:
            if buddies_to_install[buddy]:
                copytree("./{0}.py".format(buddy), "{0}/BuddySuite/{1}.py".format(install_directory, buddy))
                if install_shortcuts:
                    os.symlink("{0}/BuddySuite/{1}.py".format(install_directory, buddy),
                               "/usr/local/bin/{0}".format(buddy_to_shortcut[buddy]))

    elif user_system == 'Windows':
        return
def make_config_file():
    return