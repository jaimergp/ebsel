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