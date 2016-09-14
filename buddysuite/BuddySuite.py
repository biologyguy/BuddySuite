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
from pkg_resources import Requirement, resource_filename
from configparser import ConfigParser, NoOptionError
import os
import sys
import buddysuite
import buddysuite.buddy_resources as br
import shutil


def setup():
    valid_email_blurb = "\nProviding a valid email address is recommended if accessing public databases with " \
                        "BuddySuite.\nThe maintainers of those resources may attempt to contact you before " \
                        "blocking your IP if you are not adhering to their usage limitations.\n"

    sip_blurb = "Would you like to join our Software Improvement Program?\nAnonymized usage statistics " \
                "and crash reports will be automatically transmitted to the BuddySuite developers ([y]/n): "

    install_dir = "/".join(buddysuite.__file__.split("/")[:-2])
    os.makedirs("%s/config" % install_dir, exist_ok=True)
    if os.path.isfile("%s/config/config.ini" % install_dir):
        reader = ConfigParser()
        reader.read(config)
        email = reader.get('DEFAULT', 'email')
        diagnostics = reader.getboolean('DEFAULT', 'diagnostics')
        user_hash = reader.get('DEFAULT', 'user_hash')
        links = reader.get('DEFAULT', "links")

        print('Welcome to BuddySuite!\nPrevious installation detected...\n')
        if not diagnostics:
            diagnostics = br.ask(sip_blurb)
        else:
            diagnostics = br.ask("Would you like remain in our Software Improvement Program? ([y]/n): ")

        print(valid_email_blurb)
        email_update = input("Email address (currently '%s'): " % email)
        email = email_update if email_update not in ['', email] else email

    else:
        print('Welcome to BuddySuite!\nTo configure your installation, please answer the following questions:\n')

        diagnostics = br.ask(sip_blurb)

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


def uninstall():
    if br.ask("Are you sure you want to completely remove BuddySuite from your system y/[n]? ", default="no"):
        install_dir = "/".join(buddysuite.__file__.split("/")[:-2])

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
    import SeqBuddy
    import AlignBuddy
    import PhyloBuddy
    import DatabaseBuddy
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

    parser = argparse.ArgumentParser(prog="BuddySuite.py", formatter_class=fmt, add_help=False, usage=argparse.SUPPRESS,
                                     description='''\
\033[1mBuddySuite\033[m
  Do fun stuff with biological data files. Seriously, biological data is fun stuff :)

''')

    parser.add_argument('-setup', help='Show module version #s', action='store_true')
    parser.add_argument('-uninstall', help="List all tools", action='store_true')
    parser.add_argument('-versions', help='Show module version #s', action='store_true')
    parser.add_argument('-tools', help="List all BuddySuite tools", action='store_true')
    parser.add_argument('-count', help="Output number of tools available", action='store_true')

    in_args = parser.parse_args()

    print("%s" % "/".join(buddysuite.__file__.split("/")[:-2]))

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
