#!/usr/bin/env python
# -*- coding:utf-8 mode:python; tab-width:4; indent-tabs-mode:nil; py-indent-offset:4 -*-
##

import os
import string
import sys
from structures import BasisSetEntry

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

    def numericize(self, line, numeric_only=False, force_float=False):
        """Split a line of text by whitespace and try to convert all
        numbers to numeric values. If force_float is True, all numbers
        will be returned as floats, otherwise integers where appropriate.

        :param line: input line of text
        :type line : str
        :param numeric_only: exclude unconverted fragments if True
        :type numeric_only : bool
        :param force_float: if True, use float for all numeric values (no ints)
        :type force_float : bool
        :return: mixed list of strings and floats
        :rtype : list
        """

        converted = []
        for piece in line.split():
            try:
                v = float(piece)
                if not force_float:
                    try:
                        v = int(piece)
                    except:
                        pass
            except ValueError:
                v = piece

            if type(v) == float or numeric_only == False:
                converted.append(v)

        return converted


    def get_element_symbol(self, atomic_number):
        """Get element symbol by atomic number.

        :param atomic_nu_mber: atomic number of element, e.g. 3 for Li
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
        containing multiple elements.
        N.B.: not for ECP data!

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

    def parse_multi_g94(self, text):
        """Parse a block of Gaussian 94 format basis set data, as used by,
        Psi4, potentially containing multiple elements.
        N.B.: not for ECP data!

        :param text: a text block of basis set data for one or more elements
        :type text : str
        :return: parsed basis set data
        :rtype : list
        """

        lower = text.lower()
        spherical_or_cartesian = ""
        if "cartesian" in lower:
            spherical_or_cartesian = "cartesian"
        elif "spherical" in lower:
            spherical_or_cartesian = "spherical"


        sections = text.strip().split("****")
        parsed = [self.parse_one_g94(s) for s in sections]
        filtered = [p for p in parsed if len(p.functions) > 0]
        if spherical_or_cartesian:
            for f in filtered:
                f.spherical_or_cartesian = spherical_or_cartesian
        return filtered

    def parse_multi_from_gaussian_log_file(self, text):
        """Parse basis set data as logged by gfinput, from a log
        file.

        :param text: the contents of a log file
        :type text : str
        :return: parsed basis set data
        :rtype : list
        """

        if "gfinput" not in text.lower():
            msg = "WARNING: did not find gfinput keyword in data. This logged data may be unsuitable.\n"
            sys.stderr.write(msg)

        #Spherical or cartesian functions?
        #(5D, 7F) is spherical
        #(6D, 10F) is cartesian
        #(6D, 7F) is ???
        #(5D, 10F) is ???
        #treat cartesianness as not-sphericalness
        if "(5D, 7F)" in text:
            spherical_or_cartesian = "spherical"
        else:
            spherical_or_cartesian = "cartesian"

        #get atomic numbers from data like this:
        #---------------------------------------------------------------------
        #Center     Atomic     Atomic              Coordinates (Angstroms)
        #Number     Number      Type              X           Y           Z
        #---------------------------------------------------------------------
        #    1          6             0        0.000000    0.000000    0.000000
        #    2          1             0        0.000000    0.000000    1.083346
        #    3          1             0        1.021388    0.000000   -0.361115
        #    4          1             0       -0.510694    0.884548   -0.361115
        #    5          1             0       -0.510694   -0.884548   -0.361115
        #---------------------------------------------------------------------

        atomnos = []
        dashcount = 0

        #looking for this line (modulo white space) as marking just-before
        #basis data
        #"Number     Number      Type              X           Y           Z"
        tlist = []
        begun = False
        for line in text.split("\n"):
            pieces = line.split()
            if "Number" in pieces and "Type" in pieces and "X" in pieces and "Y" in pieces and "Z" in pieces:
                begun = True
            elif begun:
                tlist.append(line)

        tblock = "\n".join(tlist)

        for line in tblock.split("\n"):
            if "----" in line:
                dashcount += 1
                if dashcount == 2:
                    break
            else:
                numbers = self.numericize(line)
                atomnos.append(numbers[1])

        #get basis set data starting after
        # AO basis set in the form of general basis input...
        #and ending by
        # There are     N symmetry adapted...

        begin_mark = "basis set in the form of general basis input"
        end_mark = "symmetry adapted"
        after = text.split(begin_mark)[1].split(end_mark)[0]
        bs = after.rsplit("****\n", 1)[0]

        #there will be a bit of leftover junk in basis_section 0, to be cleaned
        idx = bs.find(" 1 0")
        bs = bs[idx:]

        parsed = self.parse_multi_g94(bs)
        uniques = []
        seen = set()
        #now build a list of *unique* parsed basis data in atomic order
        for number in range(1, 106):
            if number in atomnos and number not in seen:
                seen.add(number)
                j = atomnos.index(number)
                bsd = parsed[j]
                bsd.symbol = self.elements[number][0]
                atomic_number = self.get_atomic_number(bsd.symbol)
                bsd.number = atomic_number
                bsd.spherical_or_cartesian = spherical_or_cartesian
                uniques.append(bsd)

        return uniques


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
             "scale_factor" : 1.0,
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

    def parse_one_g94(self, original_text):
        """Parse a block of Gaussian 94 atomic orbital basis set data for
        one element. N.B.: not for ECP data!

        :param original_text: a text block of basis set data for one element
        :type original_text : str
        :return: parsed basis set data
        :rtype : BasisSetEntry
        """

        #get first within-asterisks block, to strip away header
        #and ignore ECP data that might be in a second block
        if "****" in original_text:
            text = original_text.split("****")[1]
        else:
            text = original_text[:]

        d = {"spherical_or_cartesian" : "spherical",
             "element_symbol" : "",
             "element_name" : "",
             "element_number" : 0,
             "scale_factor" : 1.0,
             "functions" : []}

        if "cartesian" in text.lower():
            d["spherical_or_cartesian"] = "cartesian"
        elif "spherical" in text.lower():
            d["spherical_or_cartesian"] = "spherical"

        for line in text.split("\n"):
            lower = line.lower()
            numericized = self.numericize(line)
            types = [type(n) for n in numericized]

            #need to have an alternate version where all
            # 0.7161683735d+02 etc
            #become
            #0.7161683735e+02
            #so we can see if a line would be all-numeric after replacement

            numeric_replaced = self.numericize(lower.replace('d', 'e'))
            nr_types = [type(n) for n in numeric_replaced]

            #skip comments and blank lines
            if not lower.strip() or lower[0] in ("!",):
                pass

            #this will be the element name header, like
            #Cl     0
            #can also have an integer instead of symbol with gfprint, e.g.
            #1 0
            elif types in ([str, int], [unicode, int], [int, int]):
                try:
                    element_symbol = numericized[0].title()
                    atomic_number = self.get_atomic_number(element_symbol)
                    d["element_symbol"] = element_symbol
                    d["element_number"] = atomic_number
                except AttributeError:
                    #this was a section from gfprint and we can't actually
                    #figure out the element here
                    pass

            #this will be a line of all numerical values like
            #    933.9000000              0.399612e0-02
            elif nr_types == [float] * len(types):
                d["functions"][-1][1].append(numeric_replaced)

            #could be single string, "cartesian" or "spherical" spec
            elif nr_types == [str]:
                if numericized[0].lower() not in ("spherical", "cartesian"):
                    msg = "WARNING: unknown directive or data {}\n".format(numericized)
                    sys.stderr.write(msg)

            #this will be a line heading a group of coefficients
            #S   6   1.00
            else:
                shell_type = numericized[0]
                try:
                    scale_factor = numericized[2]
                except IndexError:
                    import ipdb; ipdb.set_trace()
                d["scale_factor"] = scale_factor
                d["functions"].append((shell_type, []))

        return BasisSetEntry(d)

    def wrap_g94_to_gbs(self, basis_list, origin):
        """Ensure that data from parse_multi_from_gaussian_log_file is
        in normalized format and joined into a form suitable for .gbs
        basis files as used by Psi4.

        :param basis_list: parsed basis set data list
        :type basis_list : list
        :param origin:
        :return: .gbs-form basis set data
        :rtype : str
        """

        z = [self.format_one_g94(x, origin) for x in basis_list]
        text = self.wrap_converted_g94(z, basis_list[0].spherical_or_cartesian)
        text = text.replace("#", "!")
        text = text.replace("!BASIS", "****\n!BASIS")
        text = text.replace("****\n****\n", "****\n")
        soc_header = "{}\n".format(basis_list[0].spherical_or_cartesian)
        final = soc_header + text
        return final

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