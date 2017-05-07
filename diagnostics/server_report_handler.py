#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Created on: Aug 28 2015

"""
This program is free software in the public domain as stipulated by the Copyright Law
of the United States of America, chapter 1, subsection 105. You may modify it and/or redistribute it
without restriction.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

name: server_report_handler.py
author: Stephen R. Bond
email: steve.bond@nih.gov
institute: Computational and Statistical Genomics Branch, Division of Intramural Research,
           National Human Genome Research Institute, National Institutes of Health
           Bethesda, MD
repository: https://github.com/biologyguy/BuddySuite
© license: None, this work is public domain

Functions used by cron on the ftp server.
CAUTION!! Do not point cron directly at this unless you want a huge security hole.
In fact, DO NOT USE THIS AT ALL UNLESS YOU REALLY KNOW WHAT YOU ARE DOING!! You've been warned.
"""

import sys
import os
from datetime import date
from hashlib import md5
from subprocess import run
import re
from urllib import request
from urllib.error import URLError, HTTPError, ContentTooShortError
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import json


def sendmail(sender, recipient, subject_line, message):
    msg = MIMEMultipart()
    msg.preamble = subject_line
    msg.add_header("From", sender)
    msg.add_header("Subject", subject_line)
    msg.add_header("To", recipient)

    msg.attach(MIMEText(message))

    smtp = smtplib.SMTP('localhost')
    smtp.starttls()
    smtp.sendmail(sender, recipient, msg.as_string())
    smtp.quit()
    return

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(prog="server_report_handler", description="",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("report_folder", help="", action="store")
    parser.add_argument("-e", "--errors", help="", action="store_true")
    parser.add_argument("-u", "--usage", help="Specify location of usage file", action="store")
    parser.add_argument("-c", "--coverage", help="", action="store")

    in_args = parser.parse_args()

    root, dirs, reports = next(os.walk(in_args.report_folder))

    if in_args.errors:
        error_reports = []
        report_hashes = {}
        file_paths = []
        email_msg = None
        try:
            raw_error_data = request.urlopen("https://raw.githubusercontent.com/biologyguy/BuddySuite/error_codes/"
                                             "diagnostics/error_codes", timeout=2)
            error_string = raw_error_data.read().decode("utf-8")  # Read downloaded file
            error_string = re.sub("#.*\n", "", error_string)
            error_json = json.loads(error_string)  # Convert JSON into a data table

            resolved = []
            for key, values in error_json.items():
                report_ver, fixed_ver = values
                if fixed_ver != "":
                    resolved.append(key)

            for report in reports:
                if re.match("error", report):
                    report = "%s/%s" % (root, report)
                    file_paths.append(report)
                    with open(report, "r", encoding="utf-8") as ifile:
                        content = ifile.read()
                        _hash = md5(re.sub("#.*\n", "", content).encode()).hexdigest()

                    if _hash in resolved:
                        continue

                    if _hash not in report_hashes:
                        report_hashes[_hash] = 1
                        error_reports.append(content)
                    else:
                        report_hashes[_hash] += 1

            if len(error_reports) > 0:
                email_msg = ""
                for report in error_reports:
                    _hash = md5(report.encode()).hexdigest()
                    email_msg += "%s\n%s\n" % (_hash, report_hashes[_hash])
                    email_msg += "%s\n//\n\n" % report

        except (URLError, HTTPError, ContentTooShortError) as err:
            email_msg = "server_report_handler.py error: Unable to download error codes file"

        if email_msg:
            try:
                subject = "BuddySuite|error_reports|%s" % date.today()
                sendmail("mailer@rf-cloning.org", "buddysuite@mail.nih.gov", subject, email_msg)

                for report in file_paths:
                    os.remove(report)

            except OSError as e:
                sys.stderr("Failed to send error report:\n%s\n" % e)

        sys.exit()

    if in_args.usage:
        file_paths = []
        email_msg = ""
        for report in reports:
            if re.match("usage", report):
                report_path = "%s/%s" % (root, report)
                file_paths.append(report_path)
                with open(report_path, "r", encoding="utf-8") as ifile:
                    content = ifile.read()
                email_msg += "%s\n" % content

        if email_msg != "":
            with open(in_args.usage, "a", encoding="utf-8") as ofile:
                ofile.write(email_msg)
            try:
                subject = "BuddySuite|usage_reports|%s" % date.today()
                sendmail("mailer@rf-cloning.org", "buddysuite@gmail.com", subject, email_msg)

                for report in file_paths:
                    os.remove(report)

            except OSError as e:
                sys.stderr("Failed to send usage report:\n%s\n" % e)

        sys.exit()

    if in_args.coverage:
        if os.path.exists(".coverage"):
            run("coveralls", shell=True)
        else:
            raise FileNotFoundError("No coverage file found.")
        sys.exit()
