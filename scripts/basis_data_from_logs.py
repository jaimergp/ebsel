#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##


"""Extract reference basis set data from log files by running a series of
dummy single-atom calculations and processing the output.

Usage:
  ./scripts/basis_data_from_logs.py name-of-quantum-chemistry-program
"""

from __future__ import print_function, absolute_import
import os
import sys

workdir = "/tmp/basis_data_from_logs"

from ..ebsel import conversion

def write_one_job(basis_set_name, element):
    filename = "{}/{}__{}.com".format(workdir, basis_set_name, element)
    deck = """#n UHF/{basis} Guess=Only gfinput

 Title

0 1
{element}          0.00000        0.00000        0.00000


""".format(basis=basis_set_name, element=element)
    with open(filename, "w") as outfile:
        outfile.write(deck)

    return filename

def run_one_job(filename, qc_exe):
    for r in ["*", "(", ")", "'", " ", ","]:
        filename = filename.replace(r, "\\" + r)
    cmd = "{0} {1}".format(qc_exe, filename)
    print(cmd)
    os.system(cmd)

def extract_one_log(filename):
    try:
        with open(filename) as infile:
            data = infile.read()
    except IOError:
        print("IOERROR!")
        data = ""

    return data

def main(qc_exe):
    failures = []
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    c = conversion.Converter()

    elements = [e[0] for e in c.elements[1:106]]

    #dunning basis sets plus "calendar" variants
    dunnings = ["cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "cc-pV5Z", "cc-pV6Z"]
    calendars = []
    for prefix in ["apr", "may", "jun", "jul", "aug", "spAug", "dAug"]:
        for d in dunnings:
            bsname = "{}-{}".format(prefix, d)
            calendars.append(bsname)

    #ugbs basis set series
    ugbs = ["UGBS"]
    for n in [1, 2, 3]:
        for code in ["P", "V", "O"]:
            bsname = "UGBS{}{}".format(n, code)
            ugbs.append(bsname)

    #various basis sets with polarization and diffuse functions
    sets = [("3-21G", "", "+"), ("6-21G", "**", ""), ("4-31G", "**", ""),
            ("6-31G", "**", "++"), ("6-311G", "**", "++")]
    combinations = set()
    for base, polar, diffuse in sets:
        for j in range(3):
            for k in range(3):
                star = polar[:j]
                plus = diffuse[:k]
                name = base.replace("G", plus + "G")
                name = name.replace("G", "G" + star)
                combinations.add(name)

    combinations = sorted(list(combinations))

    basis_names = ["STO-3G", "6-31G(d')", "6-31G(d',p')", "SV", "SVP", "TZV",
                   "TZVP", "Def2SV", "Def2SVP", "Def2SVPP", "Def2TZV",
                   "Def2TZVP", "Def2TZVPP", "Def2QZV", "Def2QZVP",
                   "Def2QZVPP", "QZVP", "MidiX", "EPR-II", "EPR-III",
                   "MTSmall", "DGDZVP", "DGDZVP2", "DGTZVP", "CBSB7",
                   "6-311G(2df,p)", "6-311+G(2df,p)", "6-311++G(2df,p)",
                   "6-31G(2df,p)"]
    basis_names += combinations
    basis_names += dunnings
    basis_names += calendars
    basis_names += ugbs

    parsed = {}
    for basis in basis_names:
        parsed[basis] = []
        for element in elements:
            job_file = write_one_job(basis, element)
            logname = job_file.replace(".com", ".log")
            if not os.path.exists(logname):
                run_one_job(job_file, qc_exe)
            data = extract_one_log(logname)

            if "syntax error" in data:
                print("Syntax error in log -- misnamed basis set?")

            try:
                pbs = c.parse_multi_from_gaussian_log_file(data)[0]
                parsed[basis].append(pbs)
            except:
                pass

        if parsed[basis]:
            destination = "{}/{}.gbs".format(workdir, basis)
            origin = "{} log files".format(basis)
            for bse in parsed[basis]:
                bse.origin = origin
            combined = c.wrap_g94_to_gbs(parsed[basis])
            with open(destination, "w") as outfile:
                outfile.write(combined)
        else:
            failures.append(basis)

    os.system("rm -f ?au-*")
    if failures:
        print("FAILURES: {}".format(failures))

if __name__ == '__main__':
    qc_exe = sys.argv[1]
    main(qc_exe)