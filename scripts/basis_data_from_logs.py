#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##


"""Extract reference basis set data from log files by running a series of
dummy single-atom calculations with gfinput.

Usage:
  scripts/basis_data_from_logs.py name-of-quantum-chemistry-program
"""

import os
import glob
import sys

workdir = "/tmp/basis_data_from_logs"

from src import conversion

def write_one_job(basis_set_name, element):
    filename = "{}/{}__{}.com".format(workdir, basis_set_name, element)
    tpl = """#n UHF/{basis} Guess=Only gfinput

 Title

0 1
{element}          0.00000        0.00000        0.00000


""".format(basis=basis_set_name, element=element)
    with open(filename, "w") as outfile:
        outfile.write(tpl)

    return filename

def run_one_job(filename, qc_exe):
    cmd = "{0} {1}".format(qc_exe, filename)
    print(cmd)
    os.system(cmd)

def extract_one_log(filename):
    logname = filename.replace(".com", ".log")
    try:
        with open(logname) as infile:
            data = infile.read()
    except IOError:
        print "IOERROR!"
        data = ""

    return data

def main(qc_exe):
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    c = conversion.Converter()

    elements = [e[0] for e in c.elements[1:106]]
    basis_names = ["6-311G", "6-311+G"]
    parsed = {}
    for basis in basis_names:
        parsed[basis] = []
        for element in elements:
            job_file = write_one_job(basis, element)
            run_one_job(job_file, qc_exe)
            data = extract_one_log(job_file)
            if "Atomic number out of range" not in data:
                pbs = c.parse_multi_from_gaussian_log_file(data)[0]
                parsed[basis].append(pbs)

        if parsed[basis]:
            destination = "{}/{}.gbs".format(workdir, basis)
            origin = "{} log files".format(basis)
            combined = c.wrap_g94_to_gbs(parsed[basis], origin)
            with open(destination, "w") as outfile:
                outfile.write(combined)

    os.system("rm Gau-*")


if __name__ == '__main__':
    qc_exe = sys.argv[1]
    main(qc_exe)