#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

import sqlite3
import re
import sys
import os
import json
import time

def checkSQLite3(db_path, fmt):

    from os.path import isfile, getsize

    # Check if db file is readable
    if not os.access(db_path, os.R_OK):
        print >>sys.stderr, "Db file %s is not readable" % (db_path)
        raise IOError

    if not isfile(db_path):
        print >>sys.stderr, "Db file %s is not... a file!" % (db_path)
        raise IOError

    if getsize(db_path) < 100:  # SQLite database file header is 100 bytes
        print >>sys.stderr, "Db file %s is not a SQLite file!" % (db_path)
        raise IOError

    with open(db_path, 'rb') as fd:
        header = fd.read(100)

    if header[:16] != 'SQLite format 3\x00':
        print >>sys.stderr, "Db file %s is not in SQLiteFormat3!" % (db_path)
        raise IOError

    # Check if the file system allows I/O on sqlite3 (lustre)
    # If not, copy on /dev/shm and remove after opening
    try:
        EMSL_local(db_path, fmt).get_list_basis_available()
    except sqlite3.OperationalError:
        print >>sys.stdrerr, "I/O Error for you file system"
        print >>sys.stderr, "Try some fixe"
        new_db_path = "/dev/shm/%d.db" % (os.getpid())
        os.system("cp %s %s" % (db_path, new_db_path))
        db_path = new_db_path
    else:
        changed = False
        return db_path, changed

    # Try again to check
    try:
        EMSL_local(db_path, fmt).get_list_basis_available()
    except:
        print >>sys.stderr, "Sorry..."
        os.system("rm -f /dev/shm/%d.db" % (os.getpid()))
        raise
    else:
        print >>sys.stderr, "Working !"
        changed = True
        return db_path, changed


def cond_sql_or(table_name, l_value):

    l = []
    dmy = " OR ".join(['%s = "%s"' % (table_name, i) for i in l_value])
    if dmy:
        l.append("(%s)" % dmy)

    return l


class EMSL_local(object):

    def __init__(self, db_path, fmt="gamess-us", debug=True):
        self.db_path = db_path
        self.fmt = fmt
        self.shells = "S P D F G H I K L M".split()
        #Per-format functions to check maximum angular momentum
        self.am_checkers = {"gamess-us" : self.check_gamess_us,
                            "nwchem" : self.check_nwchem,
                            "g94" : self.check_gaussian94}

        #Per-format functions to perform extra formatting on basis set output.
        #The lambda just returns raw data unchanged.
        self.block_wrappers = {"gamess-us" : lambda x, y: x,
                               "nwchem" : self.wrap_nwchem,
                               "g94" : self.wrap_gaussian94}
        self.debug = debug

    def check_gamess_us(self, basis_blocks):
        """GAMESS-US supports only up to I basis functions. See if any
        basis blocks have higher basis functions.

        N.B.: Prior to January 2013, GAMESS-US supported only up to G basis
        functions.

        GAMESS has a special notation used in e.g. Pople basis sets: an
        "L" basis indicates S and P basis functions both with the
        same exponent. If we encounter an L function before a D function,
        it is a special-type L function and should be treated as a composite
        of S and P.

        @param basis_blocks: blocks of basis set data
        @type basis_blocks : list
        @return: (max_basis_fn, too_large)
        @rtype : tuple
        """

        shells = set(self.shells)
        greatest = 0
        d_encountered = False

        for block in basis_blocks:
            for line in block.split("\n"):
                pieces = line.split()
                try:
                    b = pieces[0]
                    if b == "D":
                        d_encountered = True

                    #If no D function has been encountered yet, an L function
                    #should be treated as a fused "SP" function, which is P for
                    #purposes of determining the greatest function AM
                    #encountered so far
                    if b == "L" and not d_encountered:
                        b = "P"
                    index = self.shells.index(b)
                    greatest = max(greatest, index)
                except (IndexError, ValueError):
                    pass

        mbf = self.shells[greatest]
        if greatest > self.shells.index("I"):
            too_large = True
        else:
            too_large = False

        return (mbf, too_large)


    def check_nwchem(self, basis_blocks):
        """NWChem supports only up to I basis functions. See if any
        basis blocks have higher basis functions.

        N.B.: This ignores ECP data.

        @param basis_blocks: blocks of basis set data
        @type basis_blocks : list
        @return: (max_basis_fn, too_large)
        @rtype : tuple
        """

        names = ["ao basis", "cd basis", "xc basis"]
        shells = set(self.shells)
        greatest = 0

        for block_json in basis_blocks:
            block_packed = json.loads(block_json)
            for name in names:
                try:
                    block = block_packed[name]
                except KeyError:
                    continue
                for line in block.split("\n")[1:]:
                    pieces = line.split()
                    try:
                        b = pieces[1]
                        index = self.shells.index(b)
                        greatest = max(greatest, index)
                    except (IndexError, ValueError):
                        pass

        mbf = self.shells[greatest]
        if greatest > self.shells.index("I"):
            too_large = True
        else:
            too_large = False

        return (mbf, too_large)

    def check_gaussian94(self, basis_blocks):
        """Gaussian and programs using its basis set format, like Psi4,
        may support arbitrarily high angular momentum basis functions.
        So only report the maximum, never declare the value too large.

        @param basis_blocks: blocks of basis set data
        @type basis_blocks : list
        @return: (max_basis_fn, too_large)
        @rtype : tuple
        """

        shells = set(self.shells)
        greatest = 0

        for block in basis_blocks:
            for line in block.split("\n")[1:]:
                pieces = line.split()
                try:
                    b = pieces[0]
                    index = self.shells.index(b)
                    greatest = max(greatest, index)
                except (IndexError, ValueError):
                    pass

        mbf = self.shells[greatest]
        too_large = False

        return (mbf, too_large)

    def wrap_gaussian94(self, blocks, basis_name):
        """Wrap up g94 blocks with **** at head and foot.

        N.B.: basis_name is currently ignored

        @param blocks: basis set data blocks
        @type blocks : list
        @param basis_name: name of the basis set
        @return: decorated basis set data blocks
        @rtype : list
        """

        nb = []
        for block in blocks:
            lines = []
            for line in block.split("\n"):
                if line.strip():
                    lines.append(line)
            nb.append("\n".join(["****"] + lines))

        nb[-1] += "\n****"
        return nb

    def spherical_or_cartesian(self, basis_name):
        """Indicate whether a basis set should be treated as using cartesian
        or pure (spherical) basis functions. A cartesian basis set uses
        6 d functions while a spherical basis set uses 5 d functions. There
        is no way to determine whether basis sets were intended for cartesian
        or spherical basis functions, or a mixture of them, other than going
        back to the original literature. One can force spherical or cartesian
        behavior for any basis set in most programs, but the choice will affect
        energies and possibly electronic convergence, due to linear
        dependencies in the basis set.

        This function simply sees if the basis set name matches one of the
        common basis sets that have cartesian functions, or not. The available
        basis sets in the EMSL Basis Set Exchange have not been thoroughly
        reviewed. Some have been cribbed from annotations in the Psi4 basis
        library.

        http://www.gaussian.com/g_tech/g_ur/m_basis_sets.htm
        

        @param basis_name: name of the basis set to test
        @type basis_name : str
        @return: "cartesian" or "spherical"
        @rtype : str
        """

        cartesians = ["3-21G", "6-21G", "4-31G", "6-31G", "6-31G*", "6-31G**",
                      "DZ (Dunning)", "DZP (Dunning)"]

        if basis_name in cartesians:
            v = "cartesian"
        else:
            v = "spherical"

        return v

    def wrap_nwchem(self, blocks, basis_name):
        """Generate NWChem basis data sections that group different
        kinds of basis set data together.

        @param blocks: fused basis set data blocks
        @type blocks : list
        @param basis_name: name of the basis set
        @type basis_name : str
        @return: basis set data sections grouped by "ao basis," "ecp," etc.
        @rtype : list
        """

        fn_type = self.spherical_or_cartesian(basis_name)
        groups = {}
        for block_json in blocks:
            block = json.loads(block_json)
            for key in block:
                try:
                    groups[key].append(block[key])
                except KeyError:
                    groups[key] = [block[key]]

        sections = []

        #Process all basis set data except ECPs.
        for btype in ["ao basis", "cd basis", "xc basis"]:
            if btype in groups:
                joined = "\n".join(groups[btype])
                s = """basis "{0}" {1}\n{2}\nEND""".format(btype, fn_type,
                                                           joined)
                sections.append(s)

        #Process ECP data if present
        ecp = groups.get("ecp")
        if ecp:
            joined = "\n".join(ecp)
            s = """ECP\n{0}\nEND""".format(joined)
            sections.append(s)

        return sections
    
    def get_list_basis_available(self, elts=[], basis=[]):
        """
        return all the basis name who contant all the elts
        """
        
        # If not elts just get the distinct name
        # Else: 1) fetch for geting all the run_id whos satisfy the condition
        #       2) Get name, description
        #       3) Parse it

        # ~#~#~#~ #
        # I n i t #
        # ~#~#~#~ #

        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()

        # ~#~#~#~#~#~ #
        # F i l t e r #
        # ~#~#~#~#~#~ #

        if basis:
            cmd_filter_basis = " ".join(cond_sql_or("name", basis, glob=True))
        else:
            cmd_filter_basis = "(1)"

        # Not Ets
        if not elts:
            cmd = """SELECT DISTINCT name, description
                     FROM basis_tab
                     WHERE {0}"""

            cmd = cmd.format(cmd_filter_basis)

        else:

            # ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~ #
            # G e t t i n g _ B a s i s I d #
            # ~#~#~#~#~#~#~#~#~#~#~#~#~#~#~ #

            str_ = """SELECT DISTINCT basis_id
                      FROM output_tab
                      WHERE elt=? AND {0}""".format(cmd_filter_basis)

            cmd = " INTERSECT ".join([str_] * len(elts)) + ";"
            c.execute(cmd, elts)

            l_basis_id = [i[0] for i in c.fetchall()]

            # ~#~#~#~#~#~#~#~#~#~#~#~#~#~ #
            # C r e a t e _ t h e _ c m d #
            # ~#~#~#~#~#~#~#~#~#~#~#~#~#~ #

            cmd_filter_basis = " ".join(cond_sql_or("basis_id", l_basis_id))
            cmd_filter_ele = " ".join(cond_sql_or("elt", elts))

            column_to_fech = "name, description"

            filter_where = " ({}) AND ({})".format(cmd_filter_ele, cmd_filter_basis)

            cmd = """SELECT DISTINCT {0}
                     FROM output_tab
                     WHERE {1}
                     ORDER BY name""".format(column_to_fech, filter_where)
        # ~#~#~#~#~ #
        # F e t c h #
        # ~#~#~#~#~ #

        c.execute(cmd)
        info = c.fetchall()
        conn.close()

        final = [i[:] for i in info]

        return final


    def get_list_element_available(self, basis_name):

        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()

        c.execute(
            "SELECT DISTINCT elt from output_tab WHERE name=:name_us COLLATE NOCASE", {
                "name_us": basis_name})

        data = c.fetchall()

        data = [str(i[0]) for i in data]

        conn.close()
        return data

    def process_raw_data(self, l_data_raw, basis_name):
        unpacked = [b[0] for b in l_data_raw]
        validator = self.am_checkers[self.fmt]
        wrapper = self.block_wrappers[self.fmt]

        self.max_am, self.am_too_large = validator(unpacked)

        if self.am_too_large and self.debug:
            msg = "WARNING: Basis set data contains angular momentum up to {0}, which is too high for {1}\n".format(self.max_am, self.fmt)
            sys.stderr.write(msg)

        transformed = wrapper(unpacked, basis_name)
        return transformed

    def get_basis(self, basis_name, elts=None):
        #  __            _
        # /__  _ _|_   _|_ ._ _  ._ _     _  _. |
        # \_| (/_ |_    |  | (_) | | |   _> (_| |
        #                                     |
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()

        if elts:
            cmd_ele = "AND " + " ".join(cond_sql_or("elt", elts))
        else:
            cmd_ele = ""

        c.execute('''SELECT DISTINCT data from output_tab
                   WHERE name="{basis_name}" COLLATE NOCASE
                   {cmd_ele}'''.format(basis_name=basis_name,
                                       cmd_ele=cmd_ele))

        l_data_raw = c.fetchall()
        conn.close()

        l_data = self.process_raw_data(l_data_raw, basis_name)
        
        return l_data

if __name__ == "__main__":

    e = EMSL_local("EMSL.db")
    l = e.get_list_basis_available()
    for i in l:
        print i

    l = e.get_list_element_available("pc-0")
    print l

    l = e.get_basis("cc-pVTZ", ["H", "He"])
    for i in l:
        print i
