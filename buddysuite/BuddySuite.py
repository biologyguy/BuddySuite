#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: BuddySuite.py
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Description: Setup and uninstall tools
"""
import argparse
from configparser import ConfigParser, NoOptionError
import os
import re
import shutil
import random
import string

try:
    import buddysuite
    import buddy_resources as br
except ImportError:
    try:
        import buddysuite
        import buddysuite.buddy_resources as br
    except AttributeError:
        from . import buddysuite
        from . import buddy_resources as br


def setup():  # ToDo: Check permissions?
    print("\033[1mWelcome to BuddySuite!\033[m\nLet's configure your installation:\n")

    install_dir, toss = os.path.split(buddysuite.__file__)
    os.makedirs("%s%sbuddy_data" % (install_dir, os.sep), exist_ok=True)
    if not os.path.isfile("{0}{1}buddy_data{1}config.ini".format(install_dir, os.sep)):
        open("{0}{1}buddy_data{1}config.ini".format(install_dir, os.sep), "w").close()

    if not os.path.isfile("{0}{1}buddy_data{1}cmd_history".format(install_dir, os.sep)):
        open("{0}{1}buddy_data{1}cmd_history".format(install_dir, os.sep), "w").close()

    if not os.path.isfile("{0}{1}buddy_data{1}buddysuite_usage.json".format(install_dir, os.sep)):
        with open("{0}{1}buddy_data{1}buddysuite_usage.json".format(install_dir, os.sep), "w") as ofile:
            ofile.write("{}")

    config = ConfigParser()
    config.read("{0}{1}buddy_data{1}config.ini".format(install_dir, os.sep))
    options = {"email": None,
               "diagnostics": None,
               "user_hash": None,
               "shortcuts": None}

    for key in options:
        try:
            if key in ['diagnostics']:
                options[key] = config.getboolean('DEFAULT', key)
            else:
                options[key] = config.get('DEFAULT', key)
        except NoOptionError:
            pass

    # Set Email
    print("\033[1mProviding a valid email address is recommended if accessing public databases with "
          "BuddySuite.\nThe maintainers of those resources may attempt to contact you before "
          "blocking your IP if you are not adhering to their usage limitations.\033[m")

    if options['email'] in [None, "buddysuite@gmail.com"]:
        email = input("Email address (optional): ")
        email = email if re.search(r"[^@]+@[^@]+\.[^@]+", email) else "buddysuite@gmail.com"
    else:
        email = input("Update email address (currently '%s'): " % options['email'])
        email = email if re.search(r"[^@]+@[^@]+\.[^@]+", email) and email not in ['', options['email']] \
            else options['email']
    options['email'] = email

    # Set up software improvement
    print("\n\033[1mBuddySuite is able to automatically send anonymized usage statistics and crash reports to the "
          "developers as part of the software improvement program.\033[m")

    question = "join" if not options['diagnostics'] else "remain in"
    options['diagnostics'] = br.ask("Would you like to %s our Software Improvement Program? [y]/n: " % question)

    # Create user hash id   'dXruTa0qkW'
    if not options["user_hash"]:
        options['user_hash'] = "".join([random.choice(string.ascii_letters + string.digits) for _ in range(10)])

    # ToDo: set up aliases or links
    # Set up shortcuts
    # print("\n\033[1mAdding shortcuts to your system PATH can make it much quicker to call the BuddySuite tools\033[m")
    options['shortcuts'] = ""

    # Write config file
    config['DEFAULT'] = options
    with open("{0}{1}buddy_data{1}config.ini".format(install_dir, os.sep), 'w') as config_file:
        config.write(config_file)

    print("""\

\033[1mSuccess! You're all set.\033[m
    Email address:     %s
    Send diagnostics:  %s

These choices can be changed at any time by re-running setup.
Enjoy the BuddySuite!
""" % (options['email'], options['diagnostics']))
    return


def uninstall():
    if br.ask("Are you sure you want to completely remove BuddySuite from your system y/[n]? ", default="no"):
        # Need to run os.path.split() twice to get to the actual install dir
        install_dir, toss = os.path.split(buddysuite.__file__)
        install_dir, toss = os.path.split(install_dir)
        # Delete all custom shortcuts
        config = br.config_values()
        for shortcut in config["shortcuts"]:
            try:
                os.remove(shortcut)
            except FileNotFoundError:
                pass

        # Delete all gateway programs
        for buddy in ['seqbuddy', 'alignbuddy', 'phylobuddy', 'databasebuddy', 'buddysuite']:
            try:
                os.remove(shutil.which(buddy))
            except FileNotFoundError:
                pass

        # Delete the main site-packages module
        try:
            shutil.rmtree(install_dir)
        except FileNotFoundError:
            pass
    return


def version():
    from buddysuite import SeqBuddy
    from buddysuite import AlignBuddy
    from buddysuite import PhyloBuddy
    from buddysuite import DatabaseBuddy
    print("SeqBuddy: %s" % SeqBuddy.VERSION.short())
    print("AlignBuddy: %s" % AlignBuddy.VERSION.short())
    print("PhyloBuddy: %s" % PhyloBuddy.VERSION.short())
    print("DatabaseBuddy: %s\n" % DatabaseBuddy.VERSION.short())


def tools():
    print("### SeqBuddy")
    for key in sorted([key for key in br.sb_flags]):
        print(key)
    print("\n### AlignBuddy")
    for key in sorted([key for key in br.alb_flags]):
        print(key)
    print("\n### PhyloBuddy")
    for key in sorted([key for key in br.pb_flags]):
        print(key)
    print("\n### DatabaseBuddy")
    for key in sorted([key for key in br.db_flags]):
        print(key)
    print()


def count():
    print("SeqBuddy: %s" % len(br.sb_flags))
    print("AlignBuddy: %s" % len(br.alb_flags))
    print("PhyloBuddy: %s" % len(br.pb_flags))
    print("DatabaseBuddy: %s" % len(br.db_flags))
    print("Total: %s\n" % sum([len(x) for x in [br.sb_flags, br.alb_flags, br.pb_flags, br.db_flags]]))


def main():
    def fmt(prog):
        return br.CustomHelpFormatter(prog)

    parser = argparse.ArgumentParser(prog="BuddySuite.py", formatter_class=fmt, usage=argparse.SUPPRESS,
                                     description='''\
\033[1mBuddySuite\033[m
  Do fun stuff with biological data files. Seriously, biological data is fun stuff :)

''')

    parser.add_argument('-setup', help='Configure BuddySuite on your system', action='store_true')
    parser.add_argument('-uninstall', help="Completely remove the system install of BuddySuite", action='store_true')
    parser.add_argument('-versions', help='Show module version #s', action='store_true')
    parser.add_argument('-tools', help="List all BuddySuite tools", action='store_true')
    parser.add_argument('-count', help="Output number of tools available", action='store_true')

    in_args = parser.parse_args()

    if in_args.versions:
        version()
    elif in_args.tools:
        tools()
    elif in_args.count:
        count()
    elif in_args.setup:
        setup()
    elif in_args.uninstall:
        uninstall()
    return

if __name__ == '__main__':
    main()
