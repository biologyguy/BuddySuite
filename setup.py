#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, Command
from pkg_resources import Requirement, resource_filename
import os
import sys
import shutil
import random
import string
from configparser import ConfigParser, NoOptionError
from tempfile import TemporaryDirectory


class CleanCommand(Command):
    """
    Custom clean command to tidy up the project root.
    http://bit.ly/2bw7xXb
    """
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    @staticmethod
    def run():
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')


class UninstallCommand(Command):
    """
    Completely remove buddysuite
    """
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    @staticmethod
    def run():
        for buddy in ["seqbuddy", "alignbuddy", "phylobuddy", "databasebuddy"]:
            buddy = shutil.which(buddy)
            if buddy:
                os.remove(buddy)
        try:
            import buddysuite
        except ImportError:
            return

        egg = "/".join(buddysuite.__file__.split("/")[:-2])
        data_dir = "%s/buddysuite_data" % "/".join(buddysuite.__file__.split("/")[:-3])
        for _next in [egg, data_dir]:
            try:
                os.remove(_next)
            except FileNotFoundError:
                pass
        return


class SetupCommand(Command):
    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    @staticmethod
    def run():
        pwd = os.getcwd()
        tmp_dir = TemporaryDirectory()
        os.mkdir("%s/buddysuite" % tmp_dir.name)
        for _file in ["AlignBuddy.py", "DatabaseBuddy.py", "PhyloBuddy.py",
                      "SeqBuddy.py", "buddy_resources.py", "__init__.py"]:
            shutil.copyfile("buddysuite/%s" % _file, "%s/buddysuite/%s" % (tmp_dir.name, _file))
        os.chdir("%s" % tmp_dir.name)

        # Need to remove the local setup version of buddysuite from Python PATH before looking for pre-installed version
        for indx, _path in enumerate(sys.path):
            if "BuddySuite" in _path:
                del sys.path[indx]
                break

        valid_email_blurb = "\nProviding a valid email address is recommended if accessing public databases with " \
                            "BuddySuite.\nThe maintainers of those resources may attempt to contact you before " \
                            "blocking your IP if you are not adhering to their usage limitations.\n"

        sip_blurb = "Would you like to join our Software Improvement Program?\nAnonymized usage statistics " \
                    "and crash reports will be automatically transmitted to the BuddySuite developers ([y]/n): "

        try:
            import buddysuite
            config = resource_filename(Requirement.parse("buddysuite"), "config/config.ini")
            reader = ConfigParser()
            reader.read(config)
            email = reader.get('DEFAULT', 'email')
            diagnostics = reader.getboolean('DEFAULT', 'diagnostics')
            user_hash = reader.get('DEFAULT', 'user_hash')

            print('Welcome to BuddySuite!\nPrevious installation detected...\n')
            if not diagnostics:
                diagnostics = ask(sip_blurb)
            else:
                diagnostics = ask("Would you like remain in our Software Improvement Program? ([y]/n): ")

            print(valid_email_blurb)
            email_update = input("Email address (currently '%s'): " % email)
            email = email_update if email_update not in ['', email] else email

        except (ImportError, NoOptionError, KeyError):
            print('Welcome to BuddySuite!\nTo configure your installation, please answer the following questions:\n')

            diagnostics = ask(sip_blurb)

            print(valid_email_blurb)
            email = input("Email address (optional): ")

            user_hash = "".join([random.choice(string.ascii_letters + string.digits) for _ in range(10)])

        writer = ConfigParser()
        writer['DEFAULT'] = {'email': email, 'diagnostics': diagnostics, 'user_hash': user_hash}
        # 'dXruTa0qkW'
        with open('config.ini', 'w') as config_file:
            writer.write(config_file)

        shutil.rmtree("buddysuite")
        try:
            os.remove("config.ini")
        except FileNotFoundError:
            pass

        for root, dirs, foiles in os.walk("./"):
            for _dir in dirs:
                shutil.rmtree("%s/%s" % (pwd, _dir), ignore_errors=True)
                shutil.move("%s/%s" % (root, _dir), "%s/%s" % (pwd, _dir))

        # Big hack to include a buddysuite_data directory in the installation directory.
        # There may be a clean way to do this, but I haven't found it.
        if 'install' in sys.argv:
            with open("bs_data.py", "w") as ofile:
                ofile.write('''
        import buddysuite
        import os

        os.makedirs("%s/buddysuite_data" % "/".join(buddysuite.__file__.split("/")[:-3]), mode=0o777, exist_ok=True)
        ''')
            from subprocess import Popen
            Popen("python bs_data.py", shell=True).wait()

        os.chdir(pwd)


def ask(input_prompt, default="yes"):
    if default == "yes":
        yes_list = ["yes", "y", '']
        no_list = ["no", "n", "abort"]
    else:
        yes_list = ["yes", "y"]
        no_list = ["no", "n", "abort", '']

    _response = input(input_prompt)
    while True:
        if _response.lower() in yes_list:
            return True
        elif _response.lower() in no_list:
            return False
        else:
            print("Response not understood. Valid options are 'yes' and 'no'.")
            _response = input(input_prompt)


PACKAGES = [
    'buddysuite',
]

DEPENDENCIES = [
    'biopython',
    'ete3',
    'six',
    'dendropy',
    'pytest',
    'scipy',
    'numpy',
    'dill',
    'suds-py3'
]

DATA_FILES = [
    ('config', ['config.ini'])
]

CLASSIFIERS = [
    'Development Status :: 3 - Alpha',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'License :: Public Domain',
    'Natural Language :: English',
    'Operating System :: MacOS',
    'Operating System :: POSIX :: Linux',
    'Programming Language :: Python :: 3 :: only',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Medical Science Apps.',
    'Topic :: Software Development :: Libraries :: Python Modules'
]

KEYWORDS = [
    'computational biology',
    'computational genetics',
    'genetics',
    'genome',
    'phylogenetics',
    'biology'
]

ENTRY_POINTS = {
    'console_scripts': [
        'alignbuddy = buddysuite.AlignBuddy:main',
        'databasebuddy = buddysuite.DatabaseBuddy:main',
        'phylobuddy = buddysuite.PhyloBuddy:main',
        'seqbuddy = buddysuite.SeqBuddy:main',
    ]
}

setup(name='buddysuite',
      version='1.1.0',
      description='BuddySuite is a collection of command line utilities written in Python for '
                  'working with biological data.',
      author='Stephen Bond',
      author_email='steve.bond@gmail.com',
      url='https://github.com/biologyguy/BuddySuite',
      packages=PACKAGES,
      install_requires=DEPENDENCIES,
      entry_points=ENTRY_POINTS,
      license='Public Domain',
      keywords=KEYWORDS,
      data_files=DATA_FILES,
      zip_safe=False,
      cmdclass={'setup': SetupCommand, 'clean': CleanCommand, 'uninstall': UninstallCommand, 'remove': UninstallCommand}
      )
