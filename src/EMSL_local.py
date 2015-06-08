#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

import conversion
import glob
import json
import os
import sqlite3
import sys

def checkSQLite3(db_path, fmt):
    # Check if db file is readable
    if not os.access(db_path, os.R_OK):
        print >>sys.stderr, "Db file %s is not readable" % (db_path)
        raise IOError

    if not os.path.isfile(db_path):
        print >>sys.stderr, "Db file %s is not... a file!" % (db_path)
        raise IOError

    if os.path.getsize(db_path) < 100:  # SQLite database file header is 100 bytes
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
        EMSL_local(db_path, fmt).get_available_basis_sets()
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
        EMSL_local(db_path, fmt).get_available_basis_sets()
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
    def __init__(self, db_path=None, fmt="gamess-us", debug=True):
        self.fmt = fmt
        if db_path is None:
            db_path = self.db_from_format(fmt)

        self.db_path = db_path

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

    def db_from_format(self, fmt):
        """Get appropriate db_path from corresponding format.

        :param fmt: format needing a db_path, e.g. "nwchem"
        :type fmt : str
        :return: path to sqlite db file
        :rtype : str
        """

        db_map = {"gamess-us" : "db/Gamess-us.db",
                  "nwchem" : "db/NWChem.db",
                  "g94" : "db/Gaussian94.db"}
        try:
            dbfile = db_map[fmt]
            db_path = os.path.dirname(os.path.dirname(__file__)) + "/" + dbfile
        except KeyError:
            msg = "Unable to find default db for format {0}\n".format(fmt)
            sys.stderr.write(msg)
            sys.exit(1)

        db_path, db_path_changed = checkSQLite3(db_path, fmt)
        return db_path

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

    def load_basis_file(self, fmt, file_name):
        """Load and parse a single supplemental basis data
        set file.

        :param file_name: name of file to load
        :type file_name : str
        :param fmt: format to load, nwchem or g94
        :return: parsed BasisSetEntry list
        :rtype : list
        """
        c = conversion.Converter()
        parser_map = {"nwchem" : c.parse_multi_nwchem,
                      "g94" : c.parse_multi_g94}

        with open(file_name) as infile:
            file_data = infile.read()
        parsefn = parser_map[fmt]
        parsed = parsefn(file_data)

        return parsed

    def get_basis_files(self, fmt):
        """Get available basis set files for supplementing basis set data
        stored in sqlite. Also warn if extraneous files are present.

        :param fmt: format to load, nwchem or g94
        :type fmt : str
        :return: basis set file paths and names
        :rtype : list
        """

        extensions = {"nwchem" : ".nwbas",
                      "g94" : ".gbs",
                      "gamess-us" : ""}

        extension = extensions[fmt]
        if not extension:
            return []

        db_root = os.path.dirname(os.path.dirname(__file__)) + "/db/"
        for ignored in ["gamess-us"]:
            pattern = db_root + ignored + "/*"
            flist = glob.glob(pattern)
            for entry in sorted(flist):
                if not entry.endswith("README.txt"):
                    msg = "WARNING: found file {} -- no files in this directory are processed\n".format(entry)
                    sys.stderr.write(msg)

        pattern = db_root + "{}/*".format(fmt)
        flist = glob.glob(pattern)
        flist.sort()
        filtered = []

        for entry in flist:
            if entry.endswith(extension):
                name = os.path.basename(entry).split(extension)[0]
                e = {"name" : name,
                     "file" : entry}
                filtered.append(e)
            elif not entry.endswith("README.txt"):
                msg = "WARNING: found unrecognized file {} -- will not be processed\n".format(entry)
                sys.stderr.write(msg)

        return filtered

    def get_available_basis_sets_fs(self, fmt, elements=[], allowed_basis_names=[]):
        """Return all the basis set names PRESENT ON THE FILE SYSTEM that
         contain the specified elements. This function looks at the
         supplementary db/nwchem/*.nwbas or db/g94/*.gbs files to find
         matching basis set data.

         If elements is empty, just get all basis set names.
         If allowed_basis_names is set, only accept results with names in
         allowed_basis_names.

        :param fmt: format to try, nwchem or g94
        :type fmt : str
        :param elements: optional element symbols to match
        :type elements : list
        :param allowed_basis_names: optional name-filter for basis set names
        :type allowed_basis_names : list
        """

        names = []
        element_set = set([e.lower() for e in elements])

        flist = self.get_basis_files(fmt)
        for entry in sorted(flist):
            t = (entry["name"], "db/" + entry["file"].split("db/", 1)[-1])
            if entry["name"] in allowed_basis_names or not allowed_basis_names:
                if not elements:
                    names.append(t)

                # if there is an element filter, need to parse actual
                # basis set data and make sure that all requested elements
                # are present
                else:
                    parsed = self.load_basis_file(fmt, entry["file"])
                    parsed_set = set([p.symbol.lower() for p in parsed])
                    if element_set.issubset(parsed_set):
                        names.append(t)

        return names
    
    def get_available_basis_sets(self, elements=[], allowed_basis_names=[]):
        """Return all the basis set names that contain the specified elements.
         If elements is empty, just get all basis set names.
         If allowed_basis_names is set, only accept results with names in
         allowed_basis_names.

        :param elements: optional element symbols to match
        :type elements : list
        :param allowed_basis_names: optional name-filter for basis set names
        :type allowed_basis_names : list
        """

        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()

        if allowed_basis_names:
            basis_filter_clause = " ".join(cond_sql_or("name", allowed_basis_names))
        else:
            basis_filter_clause = "(1)"

        if not elements:
            cmd = """SELECT DISTINCT name, description
                     FROM basis_tab
                     WHERE {0}"""

            cmd = cmd.format(basis_filter_clause)

        else:
            q = """SELECT DISTINCT basis_id
                      FROM output_tab
                      WHERE elt=? AND {0}""".format(basis_filter_clause)

            cmd = " INTERSECT ".join([q] * len(elements)) + ";"
            c.execute(cmd, elements)

            basis_ids = [i[0] for i in c.fetchall()]

            basis_filter_clause = " ".join(cond_sql_or("basis_id", basis_ids))
            elements_filter_clause = " ".join(cond_sql_or("elt", elements))

            column_to_fech = "name, description"

            filter_where = " ({}) AND ({})".format(elements_filter_clause, basis_filter_clause)

            cmd = """SELECT DISTINCT {0}
                     FROM output_tab
                     WHERE {1}
                     ORDER BY name""".format(column_to_fech, filter_where)

        c.execute(cmd)
        info = c.fetchall()
        conn.close()

        final = [i[:] for i in info]

        #look for additional basis set data from the file system
        extra = self.get_available_basis_sets_fs(self.fmt, elements=elements,
                                                 allowed_basis_names=allowed_basis_names)

        return final + extra

    def get_available_elements_fs(self, fmt, basis_name):
        """Get the available elements from a basis set that is stored
        on the file system.

        :param fmt: format to load, nwchem or g94
        :type fmt : str
        :param basis_name: name of basis set
        :type basis_name : str
        :return: element symbols
        :rtype : list
        """

        flist = self.get_basis_files(fmt)
        filtered = [x for x in flist if x["name"] == basis_name]
        if filtered:
            parsed = self.load_basis_file(fmt, filtered[0]["file"])
            elements = [p.symbol for p in parsed]
        else:
            elements = []

        return elements

    def get_available_elements(self, basis_name):
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()

        c.execute(
            "SELECT DISTINCT elt from output_tab WHERE name=:name_us COLLATE NOCASE", {
                "name_us": basis_name})

        data = c.fetchall()
        data = [str(i[0]) for i in data]
        conn.close()

        if not data:
            for fmt in ["nwchem", "g94"]:
                data = self.get_available_elements_fs(fmt, basis_name)
                if data:
                    break
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

    def convert_from_format(self, fmt, basis_name, destination_format,
                            elements=[], bypass_db=False):
        """Fetch basis set data from original specified format and return
        it in standardized converted form appropriate to destination_format.

        :param fmt: format to load, nwchem or g94
        :param basis_name: name of the basis set
        :type basis_name : str
        :param destination_format: format to convert to
        :type destination_format : str
        :param elements: elements that need basis data
        :type elements : list
        :param bypass_db: if True, look for data only on file system
        :type bypass_db : bool
        :return: basis set data for one or more elements
        :rtype : list
        """

        completed = []
        el = EMSL_local(fmt=fmt, debug=False)
        c = conversion.Converter()

        converters = {"nwchem" : c.format_one_nwchem,
                      "gamess-us" : c.format_one_gamess_us,
                      "g94" : c.format_one_g94}
        wrappers = {"nwchem" : c.wrap_converted_nwchem,
                    "gamess-us" : c.wrap_converted_gamess_us,
                    "g94" : c.wrap_converted_g94}
        dbnames = {"nwchem" : "db/NWChem.db",
                   "g94" : "db/Gaussian94.db"}
        parsers = {"nwchem" : c.parse_one_nwchem,
                   "g94" : c.parse_one_g94}

        try:
            converter = converters[destination_format]
            wrapper = wrappers[destination_format]
            parser = parsers[fmt]
            dbname = dbnames[fmt]
        except KeyError:
            raise ValueError("No defined conversion for {}".format(destination_format))

        if not bypass_db:
            if not elements:
                elements = el.get_available_elements(basis_name)

            for element in elements:
                basis = "\n".join(el.get_basis(basis_name, [element]))
                parsed = parser(basis)
                converted = converter(parsed, dbname)
                completed.append(converted)

        #either no data was found in the database or we are deliberately
        #bypassing the database to force use of basis data from file system
        if not completed:
            bs = self.get_available_basis_sets_fs(fmt, allowed_basis_names=[basis_name])

            if bs:
                flist = self.get_basis_files(fmt)
                basfile = [x["file"] for x in flist if x["name"] == basis_name][0]
                parsed = self.load_basis_file(fmt, basfile)
                if not elements:
                    elements = [p.symbol for p in parsed]

                selected = [p for p in parsed if p.symbol in elements]
                origin = bs[0][1]
                converted = [converter(p, origin) for p in selected]
                wrapped = wrapper(converted, selected[0].spherical_or_cartesian)

        else:
            wrapped = wrapper(completed, parsed.spherical_or_cartesian)

        return [wrapped]

    def get_basis(self, basis_name, elements=[], convert_from="", bypass_db=False):
        """Get basis data for named basis set. If elts is empty, all elements
         in the named basis set will be returned. If convert_from is set,
         the basis set data will first be read in the convert_from format
         and transformed to the self.fmt for output.

         If basis set data is not found in the matching native-format
         database, it will be sought in the supplemental db/nwchem/*.nwbas
         files or db/g94/*.gbs files and auto-converted if there is a match.

         The primary basis set databases can be bypassed, forcing data to load
         from the file system basis entries, if bypass_db is set to True.

        :param basis_name: name of the basis set
        :type basis_name : str
        :param elements: elements that need basis data
        :type elements : list
        :param convert_from: optional format to first convert from
        :type convert_from : str
        :param bypass_db: if True, ignore data stored in sqlite3 database
        :type bypass_db : bool
        :return: basis set data for one or more elements
        :rtype : list
        """

        processed = []
        #conversions from nwchem and g94 available presently
        if convert_from in ("nwchem", "g94"):
            return self.convert_from_format(convert_from, basis_name, self.fmt,
                                            elements=elements,
                                            bypass_db=bypass_db)
        elif convert_from:
            raise NotImplementedError("Conversion from {} not implemented".format(convert_from))

        l_data_raw = []

        if not bypass_db:
            conn = sqlite3.connect(self.db_path)
            c = conn.cursor()
            if elements:
                cmd_ele = "AND " + " ".join(cond_sql_or("elt", elements))
            else:
                cmd_ele = ""

            query = """SELECT DISTINCT data from output_tab
                       WHERE name="{basis_name}" COLLATE NOCASE
                       {cmd_ele}""".format(basis_name=basis_name,
                                           cmd_ele=cmd_ele)
            c.execute(query)

            l_data_raw = c.fetchall()
            conn.close()
            if l_data_raw:
                processed = self.process_raw_data(l_data_raw, basis_name)

        #no results from db, so try supplemenal filesystem data
        if not l_data_raw:
            #convert from nwchem by default, also if it's explicitly chosen
            if not convert_from or convert_from == "nwchem":
                processed = self.convert_from_format("nwchem", basis_name,
                                                     self.fmt,
                                                     elements=elements,
                                                     bypass_db=True)

            #convert from g94 if explicitly chosen or nwchem fallback failed
            if convert_from == "g94" or not processed:
                processed = self.convert_from_nwchem_fs("g94", basis_name, self.fmt, elements=elements)

        return processed
        
        return l_data

if __name__ == "__main__":

    e = EMSL_local("EMSL.db")
    l = e.get_available_basis_sets()
    for i in l:
        print i

    l = e.get_available_elements("pc-0")
    print l

    l = e.get_basis("cc-pVTZ", ["H", "He"])
    for i in l:
        print i
