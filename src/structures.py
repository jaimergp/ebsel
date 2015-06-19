#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

from collections import OrderedDict

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
        self.scale_factor = basis_dict["scale_factor"]
        self.basis_type = basis_dict.get("basis_type", "ao basis")
        self.origin = basis_dict.get("origin", "NO ORIGIN SUPPLIED")

    def __eq__(self, other):
        """Compare BasisSetEntries. Entries will compare as equal if
        they have the same shell structure and all the numeric data is
        equal to within 1 part per million.

        Individual shell entries to be compared look like
        ('S', [[71.61683735, 0.1543289673], [13.04509632, 0.5353281423], [3.53051216, 0.4446345422]])
        ('S', [[71.6168373, 0.154329], [13.0450963, 0.5353281], [3.5305122, 0.4446345]])

        :param other: other BSE to compare to
        :type other : BasisSetEntry
        :return: True if BSEs are equal, else False
        :rtype : bool
        """

        max_deviation = 1.0 / 1000000
        upper = 1.0 + max_deviation
        lower = 1.0 - max_deviation

        equal = True
        if repr(self) != repr(other):
            equal = False

        try:
            if len(self.functions) == len(other.functions):
                for j in range(len(self.functions)):
                    f1 = self.functions[j][1]
                    f2 = other.functions[j][1]
                    for k in range(len(f1)):
                        p1 = f1[k]
                        p2 = f2[k]
                        for m in range(len(p1)):
                            ratio = p1[m] / p2[m]
                            if ratio > upper or ratio < lower:
                                equal = False
        except IndexError:
            equal = False

        return equal

    def __ne__(self, other):
        """Inequality comparison. This is just the logical inverse of equality.
        """

        return not (self == other)


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

        #fuse together shells that share exponents
        D = PrettyOrderedDict()
        for shell, lists in ffl:
            for outer in lists:
                exps = [k[0] for k in outer]
                key = tuple([shell] + exps)
                try:
                    D[key].append(outer)
                except KeyError:
                    D[key] = [outer]

        fused = []
        for k, v in D.items():
            fused.append((k[0], v))

        return fused

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

    def format_one_nwchem(self, basis_data, origin):
        """Format one block of basis data and tag it with an
        origin comment.

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

    def format_one_gamess_us(self, basis_data, origin):
        """Format one block of basis data and tag it with an
        origin comment.

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

    def format_one_g94(self, basis_data, origin):
        """Format one block of basis data and tag it with an
        origin comment.

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
                fns.append("{}   {}   {:.1f}".format(shell, len(outer),
                                                     basis_data.scale_factor))
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

    def format_as_nwchem(self):
        return self.format_one_nwchem(self, self.origin)

    def format_as_gamess_us(self):
        return self.format_one_gamess_us(self, self.origin)

    def format_as_g94(self):
        return self.format_one_g94(self, self.origin)
