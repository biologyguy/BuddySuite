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


def setup():
    # Need to remove the local setup version of buddysuite from Python PATH before looking for pre-installed version
    for indx, _path in enumerate(sys.path):
        if "BuddySuite" in _path or "buddysuite" in _path:
            sys.path[indx] = ''
    sys.path = list(set(sys.path))

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