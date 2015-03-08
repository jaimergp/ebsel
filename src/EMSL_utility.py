# -*- coding: utf-8 -*-

import sqlite3
import re
import sys
import os
import time

debug = True


def checkSQLite3(db_path):

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
        EMSL_local(db_path=db_path).get_list_basis_available()
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
        EMSL_local(db_path=db_path).get_list_basis_available()
    except:
        print >>sys.stderr, "Sorry..."
        os.system("rm -f /dev/shm/%d.db" % (os.getpid()))
        raise
    else:
        print >>sys.stderr, "Working !"
        changed = True
        return db_path, changed


def install_with_pip(name):

    ins = False
    d = {'y': True,
         'n': False}

    while True:
        choice = raw_input('Do you want to install it ? [Y/N]')
        try:
            ins = d[choice.lower()]
            break
        except:
            print "not a valid choice"

    if ins:
        try:
            import pip
            pip.main(['install', name])
        except:
            print "You need pip, (http://pip.readthedocs.org/en/latest/installing.html)"
            sys.exit(1)


def cond_sql_or(table_name, l_value):

    l = []
    dmy = " OR ".join(['%s = "%s"' % (table_name, i) for i in l_value])
    if dmy:
        l.append("(%s)" % dmy)

    return l


class EMSL_dump:

    def __init__(self, db_path=None, format="gamess-us", contraction="True"):
        self.db_path = db_path
        self.format = format
        self.contraction = str(contraction)
        try:
            import requests
        except:
            print "You need the requests package"
            install_with_pip("requests")
        finally:
            self.requests = requests

        self.format_dict = {"g94": {"name" : "Gaussian94", "parser" : None},
                            "gamess-us": {"name" : "GAMESS-US",
                                          "parser" : self.parse_basis_data_gamess_us},
                            "gamess-uk": {"name" : "GAMESS-UK", "parser" : None},
                            "turbomole": {"name" : "Turbomole", "parser" : None},
                            "tx93": {"name" : "TX93", "parser" : None},
                            "molpro": {"name" : "Molpro", "parser" : "None"},
                            "molproint": {"name" : "MolproInt", "parser" : None},
                            "hondo": {"name" : "Hondo", "parser" : None},
                            "supermolecule": {"name" : "SuperMolecule",
                                              "parser" : None},
                            "molcas": {"name" : "Molcas", "parser" : None},
                            "hyperchem": {"name" : "HyperChem", "parser" : None},
                            "dalton": {"name" : "Dalton", "parser" : None},
                            "demon-ks": {"name" : "deMon-KS", "parser" : None},
                            "demon2k": {"name" : "deMon2k", "parser" : None},
                            "aces2": {"name" : "AcesII", "parser" : None},
                            "nwchem": {"name" : "NWChem",
                                       "parser" : self.parse_basis_data_nwchem}
                        }


    def get_list_format(self):
        """List all the format available in EMSL"""
        return self.format_dict

    def set_db_path(self, path):
        """Define the database path"""
        self.db_path = path

    def get_dict_ele(self):
        """A dict of element"""
        elt_path = os.path.dirname(sys.argv[0]) + "/src/elts_abrev.dat"
        if not os.path.exists(elt_path):
            elt_path = "src/elts_abrev.dat"

        with open(elt_path, "r") as f:
            data = f.readlines()

        dict_ele = dict()
        for i in data:
            l = i.split("-")
            dict_ele[l[1].strip().lower()] = l[2].strip().lower()
        return dict_ele

    def dwl_basis_list_raw(self):
        print "Download all the name available in EMSL. It can take some time.",
        sys.stdout.flush()

        """Download the source code of the iframe who contains the list of the basis set available"""

        url = "https://bse.pnl.gov/bse/portal/user/anon/js_peid/11535052407933/panel/Main/template/content"
        if debug:
            import cPickle as pickle
            dbcache = 'db/cache'
            if not os.path.isfile(dbcache):
                page = self.requests.get(url).text
                file = open(dbcache, 'w')
                pickle.dump(page, file)
            else:
                file = open(dbcache, 'r')
                page = pickle.load(file)
            file.close()

        else:
            page = self.requests.get(url).text

        print "Done"
        return page

    def bl_raw_to_array(self, data_raw):
        """Parse the raw html to create a basis set array whith all the info:
        url, name, description

        Explanation of tuple data from 'tup' by index:

        0  - path to xml file
        1  - basis set name
        2  - categorization: "dftcfit", "dftorb", "dftxfit", "diffuse",
             "ecporb","effective core potential", "orbital", "polarization",
             "rydberg", or "tight"
        3  - parameterized elements by symbol e.g. '[H, He, B, C, N, O, F, Ne]'
        4  - curation status; only 'published' is trustworthy
        5  - boolean: has ECP
        6  - boolean: has spin
        7  - last modified date
        8  - name of primary developer
        9  - name of contributor
        10 - human-readable summary/description of basis set
        """

        known_elements = self.get_dict_ele().keys()

        d = {}

        for line in data_raw.split('\n'):
            if "new basisSet(" in line:
                b = line.find("(")
                e = line.find(");")

                s = line[b + 1:e]

                tup = eval(s)

                #non-published (e.g. rejected) basis sets should be ignored
                if tup[4] != "published":
                    continue
                xml_path = tup[0]
                name = tup[1]

                raw_elts = re.sub('[["\ \]]', '', tup[3]).split(',')
                #filter out weird elements from exotic basis sets, e.g. Uuu
                elts = [e for e in raw_elts if e.lower() in known_elements]

                des = re.sub('\s+', ' ', tup[-1])

                if "-ecp" in xml_path.lower():
                    continue
                d[name] = [name, xml_path, des, elts]

        """Tric for the unicity of the name"""
        array = [d[key] for key in d]

        array_sort = sorted(array, key=lambda x: x[0])
        print len(array_sort), "basisset will be download"

        return array_sort

    def parse_basis_data_nwchem(self, data, name, description, elements):
        """Parse the NWChem basis data raw html to get a nice tuple.

        N.B.: Currently ignores ECP data!

        @param data: raw HTML from BSE
        @type data : unicode
        @param name: basis set name
        @type name : str
        @param des: basis set description
        @type des : str
        @param elements: element symbols e.g. ['H', 'C', 'N', 'O', 'Cl']
        @type elements : list
        @return: (name, description, data-per-element)
        @rtype : tuple
        """

        d = []

        begin_marker = """BASIS "ao basis" PRINT"""
        end_marker = "END"
        begin = data.find(begin_marker)
        end = data.find(end_marker)
        
        if begin == -1:
            if debug:
                print(data)
            raise ValueError("No basis set data found while attempting to process {0} ({1})".format(name, description))

        trimmed = data[begin+len(begin_marker) : end-len(end_marker)].strip()
        chunks = []
        lines = []

        #group lines of data delimited by #BASIS SET... into per-element chunks
        for line in trimmed.split("\n"):
            if line.startswith("#BASIS SET"):
                if lines:
                    chunks.append(lines)
                lines = [line]
            else:
                lines.append(line)

        #handle trailing chunk that is not followed by another #BASIS SET...
        if lines and lines != chunks[-1]:
            chunks.append(lines)

        #join lines back into solid text blocks
        chunks = ["\n".join(c) for c in chunks]

        #check each block for element
        unused_elements = set(elements)
        for chunk in chunks:
            #get first 3 chars of second line in block
            symbol = chunk.split("\n")[1][:3].strip()
            unused_elements.remove(symbol)

        if unused_elements:
            msg = "Warning: elements {0} left over for {1}".format(list(unused_elements), name)
            print(msg)

        return (name, description, chunks)
        
    def parse_basis_data_gamess_us(self, data, name, des, elts):
        """Parse the GAMESS-US basis data raw html to get a nice tuple"""

        d = []

        b = data.find("$DATA")
        e = data.find("$END")
        if (b == -1 or data.find("$DATA$END") != -1):
            if debug:
                print data
            raise Exception("WARNING not DATA")
        else:
            #ELEM replacements are required by buggy CRENBL ECP basis set
            replacements = [("PHOSPHOROUS", "PHOSPHORUS"),
                            ("D+", "E+"),
                            ("D-", "E-"),
                            ("ELEMTN113", "Ununtrium".upper()),
                            ("ELEMENT115", "Ununpentium".upper()),
                            ("ELEMENT117", "Ununseptium".upper())]

            for old, new in replacements:
                data = data.replace(old, new)

            split_data = data[b + 5:e - 1].split('\n\n')

            dict_ele = self.get_dict_ele()

            for (elt, data_elt) in zip(elts, split_data):

                elt_long_th = dict_ele[elt.lower()]
                elt_long_exp = data_elt.split()[0].lower()

                if "$" in data_elt:
                    if debug:
                        print "Eror",
                    raise Exception("WARNING bad split")

                if elt_long_th == elt_long_exp:
                    d.append([elt, data_elt.strip()])
                else:
                    if debug:
                        print "th", elt_long_th
                        print "exp", elt_long_exp
                        print "abv", elt
                    sys.stderr.write("WARNING not good ELEMENT\n")

        return [name, des, d]

    def create_sql(self, list_basis_array):
        """Create the sql from the list of basis available data"""

        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()

        c.execute('''CREATE TABLE basis_tab(
                            basis_id INTEGER PRIMARY KEY AUTOINCREMENT,
                                name text,
                         description text,
                                UNIQUE(name)
                  );''')

        c.execute('''CREATE TABLE data_tab(
                           basis_id INTEGER,
                                elt TEXT,
                               data TEXT,
                    FOREIGN KEY(basis_id)
                    REFERENCES basis_tab(basis_id)
                    );''')

        c.execute(''' CREATE VIEW output_tab AS
                        SELECT basis_id,
                               name,
                               description,
                               elt,
                               data
                        FROM   basis_tab
                NATURAL JOIN   data_tab
                    ''')

        conn.commit()

        import Queue
        import threading

        num_worker_threads = 1
        attemps_max = 20

        q_in = Queue.Queue(num_worker_threads)
        q_out = Queue.Queue(num_worker_threads)

        def worker():
            """get a Job from the q_in, do stuff, when finish put it in the q_out"""
            parser_method = self.format_dict[self.format]["parser"]
            if parser_method is None:
                raise NotImplementedError("No parser currently available for {0} data".format(self.format))
                
            while True:
                basis_data = []
                name, path_xml, des, elts = q_in.get()

                url = "https://bse.pnl.gov:443/bse/portal/user/anon/js_peid/11535052407933/action/portlets.BasisSetAction/template/courier_content/panel/Main/"
                url += "/eventSubmit_doDownload/true"

                params = {'bsurl': path_xml, 'bsname': name,
                          'elts': " ".join(elts),
                          'format': self.format,
                          'minimize': self.contraction}

                attemps = 0
                while attemps < attemps_max:
                    m = "URL: {0} params {1} attempt {2}".format(url, params,
                                                                 attemps)
                    print m
                    text = self.requests.get(url, params=params).text

                    #begin
                    import codecs
                    bsfname = name.replace("+", "p").replace(" ", "_").replace("*", "s").replace("/", "-") + ".html"
                    f = codecs.open(bsfname, "w", "utf-8")
                    f.write(text)
                    f.close()
                    break
                    #end
                    
                    try:
                        basis_data = parser_method(text, name, des, elts)
                    except:
                        time.sleep(0.1)
                        attemps += 1
                    else:
                        break

                """try:
                    q_out.put(basis_data)
                except:
                    if debug:
                        print "Fail on q_out.put", basis_data
                    raise
                else:
                    q_in.task_done()"""


        def enqueue():
            for [name, path_xml, des, elts] in list_basis_array:
                print "PARAMS", [name, path_xml, des, elts]
                q_in.put([name, path_xml, des, elts])

            return 0

        t = threading.Thread(target=enqueue)
        t.daemon = True
        t.start()

        for i in range(num_worker_threads):
            t = threading.Thread(target=worker)
            t.daemon = True
            t.start()

        nb_basis = len(list_basis_array)

        for i in range(nb_basis):
            name, des, d = q_out.get()
            q_out.task_done()

            try:
                c.execute(
                    "INSERT INTO basis_tab(name,description) VALUES (?,?)", [
                        name, des])
                conn.commit()
            except sqlite3.IntegrityError:
                print '{:>3}'.format(i + 1), "/", nb_basis, name, "fail"

            id_ = c.lastrowid
            try:
                c.executemany(
                    "INSERT INTO data_tab VALUES (?,?,?)", [
                        [id_] + k for k in d])
                conn.commit()

                print '{:>3}'.format(i + 1), "/", nb_basis, name

            except:
                print '{:>3}'.format(i + 1), "/", nb_basis, name, "fail"
                raise

        conn.close()

        q_in.join()

    def new_db(self):
        """Create new_db from scratch"""

        _data = self.dwl_basis_list_raw()
        array_basis = self.bl_raw_to_array(_data)
        del _data

        self.create_sql(array_basis)


class EMSL_local:

    def __init__(self, db_path=None):
        self.db_path = db_path

    def get_list_basis_available(self, elts=[]):

        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()

        if not elts:

            c.execute("""SELECT DISTINCT name,description, LENGTH(data)- LENGTH(REPLACE(data, X'0A', ''))
                         FROM output_tab""")
            data = c.fetchall()

        else:
            cmd = ["""SELECT name,description, LENGTH(data)- LENGTH(REPLACE(data, X'0A', ''))
                      FROM output_tab WHERE elt=?"""] * len(elts)
            cmd = " INTERSECT ".join(cmd) + ";"

            c.execute(cmd, elts)
            data = c.fetchall()

        data = [i[:] for i in data]

        conn.close()

        return data

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

    def get_basis(self, basis_name, elts=None, with_l=False):

        import re

        def get_list_type(l_line):
            l = []
            for i, line in enumerate(l_line):

                m = re.search(p, line)
                if m:
                    l.append([m.group(1), i])
                    try:
                        l[-2].append(i)
                    except IndexError:
                        pass

            l[-1].append(i + 1)
            return l

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

        # |_|  _. ._   _| |  _    || | ||
        # | | (_| | | (_| | (/_      |_
        #

        p = re.compile(ur'^(\w)\s+\d+\b')

        l_data = []

        for data_raw in l_data_raw:

            basis = data_raw[0].strip()

            l_line_raw = basis.split("\n")

            l_line = [l_line_raw[0]]

            for symmetry, begin, end in get_list_type(l_line_raw):

                if not(with_l) and symmetry in "L":

                    body_s = []
                    body_p = []

                    for i_l in l_line_raw[begin + 1:end]:

                        a = i_l.split()

                        common = "{:>3}".format(a[0])
                        common += "{:>15.7f}".format(float(a[1]))

                        tail_s = common + "{:>23.7f}".format(float(a[2]))
                        body_s.append(tail_s)

                        tail_p = common + "{:>23.7f}".format(float(a[3]))
                        body_p.append(tail_p)

                    l_line += [l_line_raw[begin].replace("L", "S")]
                    l_line += body_s

                    l_line += [l_line_raw[begin].replace("L", "P")]
                    l_line += body_p
                else:
                    l_line += l_line_raw[begin:end]

            l_data.append("\n".join(l_line))

        return l_data

if __name__ == "__main__":

    e = EMSL_local(db_path="EMSL.db")
    l = e.get_list_basis_available()
    for i in l:
        print i

    l = e.get_list_element_available("pc-0")
    print l

    l = e.get_basis("cc-pVTZ", ["H", "He"])
    for i in l:
        print i
