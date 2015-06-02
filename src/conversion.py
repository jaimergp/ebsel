#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

from collections import OrderedDict
import os
import string

class PrettyOrderedDict(OrderedDict):
    def __str__(self):
        tpl = "{} : {}, " * len(self)
        tpl = tpl[:-2]

        tokens = []
        for key, value in self.items():
            tokens.append(repr(key))
            tokens.append(repr(value))

        formatted = tpl.format(*tokens)
        return "{" + formatted + "}"

class BasisSetEntry(object):
        def __init__(self, basis_dict):
            self.symbol = basis_dict["element_symbol"]
            self.number = basis_dict["element_number"]
            self.name = basis_dict["element_name"]
            self.spherical_or_cartesian = basis_dict["spherical_or_cartesian"]
            self.functions = basis_dict["functions"]

        def reformat_functions(self):
            return self._reformat_functions(self.functions)

        def _reformat_functions(self, function_list):
            """These are equivalent:

            c1 e1 e2
            c2 e1 e2

            c1 e1
            c2 e1
            c1 e2
            c2 e2

            Reformat nested function lists of the first layout to produce
            the second layout.

            e.g. take each
            ('S', [[192.1714, 0.5289 0.2731],
                   [86.1207, 0.1208, 0.0303]])

            and make it
            ('S', [[192.1714, 0.5289],
                   [86.12, 0.1208],
                   [192.1714, 0.2731],
                   [86.1207, 0.0303]])

            DO NOT do this for SP functions ("L" functions in GAMESS
            terminology). The programs demand a 3 column format for
            SP functions.

            :param function_list: functions to reformat
            :type function_list : list
            :return: restructured function list
            :rtype : list
            """

            ffl = []

            for shell, lists in function_list:
                if shell == "SP":
                    ffl.append((shell, [lists]))
                else:
                    n_columns = len(lists[0]) - 1
                    columns = [list() for c in range(n_columns)]

                    for fn in lists:
                        for j, e in enumerate(fn[1:]):
                            entry = [fn[0], e]
                            columns[j].append(entry)

                    ffl.append((shell, columns))

            return ffl

        @property
        def functions_per_shell(self):
            """Return the number of basis functions by shell name.

            e.g. for cc-pVDZ Li the return value will be
            {"s" : 3, "p" : 2, "d" : 1}

            :return: count of basis functions for each shell name
            :rtype : dict
            """
            try:
                return self._functions_per_shell
            except AttributeError:
                pass

            d = {}
            for shell, values in self.functions:
                n = len(values[0]) - 1
                try:
                    d[shell] += n
                except KeyError:
                    d[shell] = n

            #order the shells by angular momentum instead of alphabetically
            #for display
            shells = ["S", "P", "SP", "D", "F", "G", "H", "I", "K", "L", "M"]
            D = PrettyOrderedDict()
            for s in shells:
                if s in d:
                    D[s] = d[s]

            self._functions_per_shell = D
            return D

        def __repr__(self):
            s = "<{0} {1} {2}>".format(self.symbol,
                                       self.spherical_or_cartesian,
                                       self.functions_per_shell)
            return s

class Converter(object):
    def __init__(self):
        self._prepare_element_data()

    def _prepare_element_data(self):
        """Prepare elements data in atomic number order, with symbol and
        name for each element.
        """

        datfile = os.path.dirname(__file__) + "/elts_abrev.dat"
        if not os.path.exists(datfile):
            datfile = "src/elts_abrev.dat"

        with open(datfile, "r") as f:
            data = f.readlines()

        #fill zeroth entry with dummy tuple so we can index elements directly
        #by atomic number
        self.elements = [("", "")]

        for i in data:
            l = i.split("-")
            symbol = l[1].strip()
            element = l[2].strip().upper()
            self.elements.append((symbol, element))

    def get_element_symbol(self, atomic_number):
        """Get element symbol by atomic number.

        :param atomic_number: atomic number of element, e.g. 3 for Li
        :type atomic_number : int
        :return: element's symbol
        :rtype : str
        """

        return self.elements[atomic_number][0]

    def get_element_name(self, atomic_number):
        """Get element name by atomic number.

        :param atomic_number: atomic number of element, e.g. 3 for LITHIUM
        :type atomic_number : int
        :return: element's name
        :rtype : str
        """

        return self.elements[atomic_number][1]

    def get_atomic_number(self, symbol):
        """Get element atomic number by symbol.

        :param symbol: symbol of element, e.g. Li for 3
        :type symbol : str
        :return: element's atomic number
        :rtype : int
        """

        slower = [j[0].lower() for j in self.elements]
        i = slower.index(symbol.lower())
        return i

    def parse_multi_nwchem(self, text):
        """Parse a block of NWChem atomic orbital basis set data potentially
        containing multiple elements. N.B.: not for ECP data!

        :param text: a text block of basis set data for one or more elements
        :type text : str
        :return: parsed basis set data
        :rtype : list
        """

        chunks = []
        for line in text.split("\n"):
            lstrip = line.lower().strip()
            if not lstrip:
                continue
            if lstrip.startswith("basis"):
                chunks.append([line])
            else:
                #skip any comment lines encountered before first element
                if not chunks:
                    pass
                else:
                    chunks[-1].append(line)

        rejoined = ["\n".join(c) for c in chunks]
        parsed = [self.parse_one_nwchem(r) for r in rejoined]
        return parsed

    def parse_one_nwchem(self, text):
        """Parse a block of NWChem atomic orbital basis set data for
        one element. N.B.: not for ECP data!

        :param text: a text block of basis set data for one element
        :type text : str
        :return: parsed basis set data
        :rtype : BasisSetEntry
        """

        d = {"spherical_or_cartesian" : "spherical",
             "element_symbol" : "",
             "element_name" : "",
             "element_number" : 0,
             "functions" : []}

        for line in text.split("\n"):
            lower = line.lower()
            if lower.startswith("#") or not lower:
                pass

            elif lower.startswith("end") or lower.startswith("ecp"):
                break

            #this will be a line starting a basis set like
            #basis "ao basis" spherical
            elif lower.startswith("basis"):
                if "spherical" in lower:
                    d["spherical_or_cartesian"] = "spherical"
                elif "cartesian" in lower:
                    d["spherical_or_cartesian"] = "cartesian"

            #this will be a line of numerical values like
            #    508.4400000              0.4365900
            elif lower[0] in string.whitespace:
                #Python doesn't understand notation like
                #0.4137d-06
                #but
                #0.4137e-06
                #works
                lower = lower.replace("d", "e")
                values = [float(j) for j in lower.split()]
                d["functions"][-1][1].append(values)

            #this will be a line heading a group of coefficients
            #Na   SP
            else:
                symbol, shell_type = line.split()
                shell_type = shell_type.upper()
                atomic_number = self.get_atomic_number(symbol)
                element_symbol = self.get_element_symbol(atomic_number)
                d["element_symbol"] = element_symbol
                d["element_number"] = atomic_number
                d["element_name"] = self.get_element_name(atomic_number)
                d["functions"].append((shell_type, []))

        return BasisSetEntry(d)

    def format_one_nwchem(self, basis_name, basis_data, origin):
        """Format one block of basis data and tag it with an
        origin comment.

        :param basis_name: name of the basis set that provided the data
        :type basis_name : str
        :param basis_data: a standard "tall" basis set entry
        :type basis_data : BasisSetEntry
        :param origin: where the data originally came from
        :type origin : str
        :return: a formatted basis data block for NWChem
        :rtype : str
        """

        fns = []
        fps = basis_data.functions_per_shell
        contracted = []
        for key, value in fps.items():
            entry = "{}{}".format(value, key.lower())
            contracted.append(entry)

        reformatted = basis_data.reformat_functions()
        for shell, functions in reformatted:
            for outer in functions:
                fns.append("{}     {}".format(basis_data.symbol, shell))
                for vals in outer:
                    col1 = "{:.7f}".format(vals[0])
                    col2 = "{:.7f}".format(vals[1])
                    try:
                        col3 = "{:.7f}".format(vals[2])
                        pad3 = " " * (16 - col3.index("."))
                    except IndexError:
                        col3 = ""
                        pad3 = ""
                    pad1 = " " * (8 - col1.index("."))
                    pad2 = " " * (16 - col2.index("."))

                    row = pad1 + col1 + pad2 + col2 + pad3 + col3
                    fns.append(row)

        c2 = "#BASIS SET reformatted: [{}]".format(",".join(contracted))
        c3 = "#origin: {}".format(origin)

        block = "\n".join([c2, c3] + fns)
        return block

    def wrap_converted_nwchem(self, basis_set_entries, spherical_or_cartesian):
        """Wrap a list of converted basis set entries into a basis
        set data section suitable for embedding in NWChem input decks.

        :param basis_set_entries: one or more BasisSetEntry values to wrap
        :type basis_set_entries : list
        :param spherical_or_cartesian: pure or cartesian form functions
        :type spherical_or_cartesian : str
        :return: formatted basis set data section
        :rtype : str
        """

        if not basis_set_entries:
            formatted = ""
        else:
            head = """basis "ao basis" {} """.format(spherical_or_cartesian)
            formatted = "\n".join([head] + basis_set_entries + ["END"])

        return formatted

    def format_one_gamess_us(self, basis_name, basis_data, origin):
        """Format one block of basis data and tag it with an
        origin comment.

        :param basis_name: name of the basis set that provided the data
        :type basis_name : str
        :param basis_data: a standard "tall" basis set entry
        :type basis_data : BasisSetEntry
        :param origin: where the data originally came from
        :type origin : str
        :return: a formatted basis data block for GAMESS-US
        :rtype : str
        """

        fns = []
        fps = basis_data.functions_per_shell
        contracted = []
        for key, value in fps.items():
            entry = "{}{}".format(value, key.lower())
            contracted.append(entry)

        reformatted = basis_data.reformat_functions()
        for shell, functions in reformatted:
            #GAMESS calls SP shells L shells
            if shell == "SP":
                shell = "L"
            for outer in functions:
                fns.append("{}   {}".format(shell, len(outer)))
                for j, vals in enumerate(outer):
                    col0 = str(j + 1)
                    pad0 = " " * (3 - len(col0))
                    col1 = "{:.7f}".format(vals[0])
                    col2 = "{:.7f}".format(vals[1])
                    try:
                        col3 = "{:.7f}".format(vals[2])
                        pad3 = " " * (15 - col3.index("."))
                    except IndexError:
                        col3 = ""
                        pad3 = ""
                    pad1 = " " * (7 - col1.index("."))
                    pad2 = " " * (15 - col2.index("."))
                    row = pad0 + col0 + pad1 + col1 + pad2 + col2 + pad3 + col3
                    fns.append(row)

        c1 = basis_data.name.upper()
        c2 = "!BASIS SET reformatted: [{}]".format(",".join(contracted))
        c3 = "!origin: {}".format(origin)

        block = "\n".join([c1, c2, c3] + fns)
        return block

    def wrap_converted_gamess_us(self, basis_set_entries, spherical_or_cartesian):
        """Wrap a list of converted basis set entries into a basis
        set data section suitable for embedding in GAMESS-US input decks.

        :param basis_set_entries: one or more BasisSetEntry values to wrap
        :type basis_set_entries : list
        :param spherical_or_cartesian: pure or cartesian form functions
        :type spherical_or_cartesian : str
        :return: formatted basis set data section
        :rtype : str
        """

        formatted = "\n".join(basis_set_entries)
        return formatted

    def format_one_g94(self, basis_name, basis_data, origin):
        """Format one block of basis data and tag it with an
        origin comment.

        :param basis_name: name of the basis set that provided the data
        :type basis_name : str
        :param basis_data: a standard "tall" basis set entry
        :type basis_data : BasisSetEntry
        :param origin: where the data originally came from
        :type origin : str
        :return: a formatted basis data block for Gaussian 94 or compatible
        :rtype : str
        """

        fns = []
        fps = basis_data.functions_per_shell
        contracted = []
        for key, value in fps.items():
            entry = "{}{}".format(value, key.lower())
            contracted.append(entry)

        fns.append("{}     0".format(basis_data.symbol))

        reformatted = basis_data.reformat_functions()
        for shell, functions in reformatted:
            for outer in functions:
                fns.append("{}   {}   1.0".format(shell, len(outer)))
                for j, vals in enumerate(outer):
                    col1 = "{:.7f}".format(vals[0])
                    col2 = "{:.7f}".format(vals[1])
                    try:
                        col3 = "{:.7f}".format(vals[2])
                        pad3 = " " * (15 - col3.index("."))
                    except IndexError:
                        col3 = ""
                        pad3 = ""
                    pad1 = " " * (7 - col1.index("."))
                    pad2 = " " * (15 - col2.index("."))
                    row = pad1 + col1 + pad2 + col2 + pad3 + col3
                    fns.append(row)

        c2 = "#BASIS SET reformatted: [{}]".format(",".join(contracted))
        c3 = "#origin: {}".format(origin)

        block = "\n".join([c2, c3] + fns)
        return block

    def wrap_converted_g94(self, basis_set_entries, spherical_or_cartesian):
        """Wrap a list of converted basis set entries into a basis
        set data section suitable for using in Gaussian 94 or compatible
        format.

        :param basis_set_entries: one or more BasisSetEntry values to wrap
        :type basis_set_entries : list
        :param spherical_or_cartesian: pure or cartesian form functions
        :type spherical_or_cartesian : str
        :return: formatted basis set data section
        :rtype : str
        """

        formatted = "\n".join(["****"] + basis_set_entries + ["****"])
        return formatted