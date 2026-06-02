import numpy as np
import matplotlib.pyplot as plt 
import os
import math
import functools
import matplotlib as mpl
from matplotlib import cm
import xlsxwriter
import sys
from pymatgen.analysis.reaction_calculator import ComputedReaction
from pymatgen.core import Composition, Element
from pymatgen.ext.matproj import MPRester
import copy
import matplotlib.patches as mpatches

def seperate_string_number(string):
    previous_character = string[0]
    groups = []
    ##if not previous_character.isdigit(): groups.append("1")
    newword = string[0]
    par_start = 0
    par_end = 0
    #par_multiply = False
    #print(f"compound: {string}")
    if (len(string) == 1): groups.append(string)
    else:
        for x, i in enumerate(string[1:]):
            if (i.isalpha() and i.islower()) and (previous_character.isalpha() and previous_character.isupper()):
                newword += i
            elif (i.isalpha() and i.islower()) and (previous_character.isalpha() and previous_character.islower()):
                newword += i
            elif (i.isalpha() and i.isupper()) and (previous_character.isalpha() and previous_character.islower()):
                groups.append(newword)
                groups.append("1")
                newword = i
            elif (i.isalpha() and i.isupper()) and (previous_character.isalpha() and previous_character.isupper()):
                groups.append(newword)
                groups.append("1")
                newword = i
            elif i.isnumeric() and previous_character.isnumeric():
                newword += i

            elif i.isnumeric() and (previous_character == ")"):
                for j in np.arange(par_start+1, par_end+1):
                    if groups[j].isnumeric():
                        # print(f"parentheses multiplier: j: {j}  groups[j]: {groups[j]}  i: {int(i)}")
                        groups[j] = str( int(groups[j]) * int(i) )
                        # print(f"post groups[j]: {groups[j]}")

            elif (i == "(") and previous_character.isnumeric():
                groups.append(newword)
                par_start = len(groups)
            elif (i == "(") and previous_character.isalpha():
                groups.append(newword)
                groups.append("1")
                par_start = len(groups)
            elif (i == ")"):
                par_end = len(groups)
                groups.append(newword)
            elif (previous_character == "("):
                newword = i
            else:
                groups.append(newword)
                newword = i

            previous_character = i

            if x == len(string) - 2:
                groups.append(newword)
                newword = ''

    #print(f"groups: {groups}\n")
    if not groups[-1].isdigit(): groups.append("1")
    
    #print(f"groups: {groups}")

    return groups


def get_most_stable_entry(formula, all_compounds):
    a = MPRester("GV69Do38gPOc4uk2Y4ZjFk9Wwg9p5xgw")
    all_entries = a.get_entries_in_chemsys(all_compounds)
    relevant_entries = [entry for entry in all_entries if entry.composition.reduced_formula == Composition(formula).reduced_formula]
    relevant_entries = sorted(relevant_entries, key=lambda e: e.energy_per_atom)
    #print(f"relevant_entries: {relevant_entries}")
    return relevant_entries[0]


def reduce_compound(compound):
    compound_split = seperate_string_number(compound)
    coefficients = np.asarray([int(num) for num in compound_split if num.isdigit()])
    elements = np.asarray([str(elem) for elem in compound_split if not elem.isdigit()])
    gcd = 0
    #print(f"compound: {compound}")
    #print(f"coefficients: {coefficients}")
    #print(f"elements: {elements}")

    while (gcd != 1):
        gcd = functools.reduce(lambda x,y:math.gcd(x,y),coefficients)
        for i in range(len(coefficients)):  coefficients[i] = int(coefficients[i] / gcd)

    coefficients_sorted = [x for _, x in sorted(zip(elements, coefficients))]
    elements = np.sort(elements)

    new_compound = ""
    for i in range(len(elements)):
        new_compound += elements[i]
        new_compound += str(coefficients_sorted[i])

    #print(f"new_compound: {new_compound}")
    #print(f"coefficients: {coefficients_sorted}")

    return new_compound


def reduce_compound_and_coeff(compound):
    compound_split = seperate_string_number(compound)
    coefficients = np.asarray([int(num) for num in compound_split if num.isdigit()])
    elements = np.asarray([str(elem) for elem in compound_split if not elem.isdigit()])
    gcd = 0
    #print(f"compound: {compound}")
    #print(f"coefficients: {coefficients}")
    #print(f"elements: {elements}")

    while (gcd != 1):
        gcd = functools.reduce(lambda x,y:math.gcd(x,y),coefficients)
        for i in range(len(coefficients)):  coefficients[i] = int(coefficients[i] / gcd)

    coefficients_sorted = [x for _, x in sorted(zip(elements, coefficients))]
    elements = np.sort(elements)

    new_compound = ""
    for i in range(len(elements)):
        new_compound += elements[i]
        new_compound += str(coefficients_sorted[i])

    #print(f"new_compound: {new_compound}")
    #print(f"coefficients: {coefficients_sorted}")

    return new_compound, gcd


def parse_reaction(react_line):
    reactant_dict = {}
    product_dict = {}

    if (react_line.split()[-1] == "eV/atom"):
        #print(f"react_line: {react_line}")
        ### initialize dictionaries for reactants and products ###
        
        ### get energy ###
        #print(f"react_line pre: {react_line}")
        react_energy = float(react_line.split()[-2])
        react_line = " ".join(react_line.split()[:-2])

    else: react_energy = None

    ### get reaction and product side ###
    react_line_split = react_line.split("->")
    reactant_side = react_line_split[0]
    product_side = react_line_split[1]

    ### split both sides along (+) to get compounds & coefficients ###
    reactant_side_comp_coeff = reactant_side.split(" + ")
    product_side_comp_coeff = product_side.split(" + ")

    for reactant in reactant_side_comp_coeff:
        reactant_split = reactant.strip().split()

        if (len(reactant_split) == 2):
            coeff = float(reactant_split[0])
            compound = reactant_split[1]

        else: 
            coeff = 1.0 # assume default value of 1 if no coefficient
            compound = reactant_split[0]
        
        #print(f"compound: {compound}")
        compound_split = seperate_string_number(compound)
        coefficients = np.asarray([int(num) for num in compound_split if num.isdigit()])
        elements = np.asarray([str(elem) for elem in compound_split if not elem.isdigit()])
  
        coefficients_sorted = [x for _, x in sorted(zip(elements, coefficients))]
        elements = np.sort(elements)

        new_compound = ""
        for i in range(len(elements)):
            new_compound += elements[i]
            new_compound += str(coefficients_sorted[i])
            #print(f"elements[i]: {elements[i]}   coefficients_sorted[i]: {coefficients_sorted[i]}")

        #print(f"reactant compound: {compound}  coeff: {coeff}")
        if new_compound not in reactant_dict: reactant_dict[new_compound] = coeff
        else:  reactant_dict[new_compound] += coeff

    for product in product_side_comp_coeff:
        product_split = product.strip().split()

        if (len(product_split) == 2):
            coeff = float(product_split[0])
            compound = product_split[1]

        else: 
            coeff = 1 # assume default value of 1 if no coefficient
            compound = product_split[0]

        #print(f"coeff: {coeff}")
        #print(f"compound: {compound}")
        compound_split = seperate_string_number(compound)
        coefficients = np.asarray([int(num) for num in compound_split if num.isdigit()])
        elements = np.asarray([str(elem) for elem in compound_split if not elem.isdigit()])
  
        coefficients_sorted = [x for _, x in sorted(zip(elements, coefficients))]
        elements = np.sort(elements)

        new_compound = ""
        for i in range(len(elements)):
            new_compound += elements[i]
            new_compound += str(coefficients_sorted[i])

        
        #print(f"product compound: {compound}  coeff: {coeff}")
        if new_compound not in product_dict: product_dict[new_compound] = coeff
        else:  product_dict[new_compound] += coeff

    return reactant_dict, product_dict, react_energy


def check_if_redox(react_dict, prod_dict):

    for reactant in react_dict:
        elem = element(reactant)
        charges = []
        for ir in elem.ionic_radii:
            charges.append(ir.charge)


def make_reaction_subscripted(react_line):
    
    ### get reaction and product side ###
    react_line_split = react_line.split("->")
    reactant_side = react_line_split[0]
    product_side = react_line_split[1]

    ### split both sides along (+) to get compounds & coefficients ###
    reactant_side_comp_coeff = reactant_side.split(" + ")
    product_side_comp_coeff = product_side.split(" + ")

    new_reactant_str = ""
    for reactant in reactant_side_comp_coeff:
        reactant_split = reactant.strip().split()

        if (len(reactant_split) == 2):
            coeff = float(reactant_split[0])
            compound = reactant_split[1]

        else: 
            coeff = 1.0 # assume default value of 1 if no coefficient
            compound = reactant_split[0]
        
        SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")  
        subscripted_compound = compound.translate(SUB)
        
        if (coeff != 1): new_reactant_str += str(coeff) + " " + subscripted_compound
        else: new_reactant_str += subscripted_compound
        
        if (reactant == reactant_side_comp_coeff[-1]): new_reactant_str += " \u2192 "
        else:  new_reactant_str += " + "

    for product in product_side_comp_coeff:
        product_split = product.strip().split()

        if (len(product_split) == 2):
            coeff = product_split[0]
            compound = product_split[1]

        else: 
            coeff = 1 # assume default value of 1 if no coefficient
            compound = product_split[0]
        
        SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")  
        subscripted_compound = compound.translate(SUB)

        if (coeff != 1): new_reactant_str += str(coeff) + " " + subscripted_compound
        else: new_reactant_str += subscripted_compound
        
        if (product == product_side_comp_coeff[-1]): pass
        else: new_reactant_str += " + "

    return new_reactant_str


def parse_volfile(vol_lines):
    ### initialize dictionaries for reactants and products ###
    compound_dict = {}

    ### split both line along (+) to get compounds & coefficients ###
    for vol_line in vol_lines:
        if ("--------" in vol_line) or ("Entries with" in vol_line) or ("Whole compound" in vol_line) or vol_line.isspace(): pass
        else:
            vol_line_split = vol_line.strip().split()
            parse_i = 0

            while (")" not in vol_line_split[parse_i]): parse_i += 1
            compound_unreduced = "".join(vol_line_split[:(parse_i+1)])
            compound_unreduced = compound_unreduced.strip(")").strip("(")
            compound = reduce_compound(compound_unreduced)
            
            if (vol_line_split[-1] == "[cm^3/mol]"): molar_mass = float(vol_line_split[-2])
            else: molar_mass = float(vol_line_split[-3])
            
            compound_dict[compound] = molar_mass

    compound_dict['Li1'] = 13.092
    compound_dict['Na1'] = 23.948

    return compound_dict


def assign_vol_changes(volfilepath, reactfilepath, reactfolder, open_vol_included):
    
    volfile = open("stable_products_volchange_new/" + volfilepath, "r")
    
    volfile_lines = [line for line in volfile]
    
    reactfile = open(f"{reactfolder}/{reactfilepath}", "r")
    reactfile_lines = [line for line in reactfile]

    print(f"reactfile: {reactfilepath}")
    
    reactfile_split = (reactfilepath.split("_"))
    elems = reactfile_split[:-2]
    elems_joined =  "_".join(elems)
    voltage = reactfile_split[-2]
    vol_pot_name = "volumes_by_voltage/" + elems_joined
    vol_pot_name += "_volumes_by_voltage.txt"
    vol_by_pot_file = open(vol_pot_name, 'a')
    
    compound_dict = parse_volfile(volfile_lines)
    print(f"compound_dict: {compound_dict}")
    vol_by_pot_file.write(f"elements: {elems_joined} at {voltage}V\n") 
    vol_by_pot_file.write(f"--------------------------------------\n") 

    reac_vol_dict = {}
    reac_energy_dict = {}
    reac_electrolyte_mol_dict = {}
    reac_strain_dict = {}
    react_open_el_ratio_dict = {}
    reac_vol_ratio_dict = {}
    
    for react_line in reactfile_lines:
        print(f"react_line: {react_line}")
        if ("--------" in react_line) or ("decomposition reactions" in react_line) or react_line.isspace(): pass
        else:
            #print(f"react_line: {react_line}")
            react_dict, prod_dict, react_energy = parse_reaction(react_line)
            print(f"react_dict: {react_dict}")
            print(f"react_energy: {react_energy}")

            ### generating a list of all the elements in reaction
            all_elems = []

            for compound, coeff in react_dict.items():
                compound_split = seperate_string_number(compound)
                elements = np.asarray([str(elem) for elem in compound_split if not elem.isdigit()])
                for elem in elements: all_elems.append(elem)

            for compound, coeff in prod_dict.items():
                compound_split = seperate_string_number(compound)
                elements = np.asarray([str(elem) for elem in compound_split if not elem.isdigit()])
                for elem in elements: all_elems.append(elem)

            all_elems = list(set(all_elems))


            reactant_vol = 0
            electrolyte_coeff = 0    
            open_el_coeff = 0
                
            print(f"compound_dict: {compound_dict}")     
            for compound, coeff in react_dict.items():
                print(f"reactant compound: {compound}")
                reduced_compound, gcd = reduce_compound_and_coeff(compound)
                if (reduced_compound not in compound_dict): 
                    KeyError("Key 'compound' is not in dictionary 'compound_dict'")
                    exit()
                elif (((compound == "Li1") or (compound == "Na1")) and (not open_vol_included)): 
                     print(f"reactant exluded: {reduced_compound}")
                     pass # include this if statement in case of omitting volume of open element
                else: 
                    molar_volume = float(gcd) * compound_dict[reduced_compound]
                    reactant_vol += float(molar_volume) * float(coeff)

                print(f"volume: {molar_volume}")  
                compound_split = seperate_string_number(compound)
                compound_elements = [str(elem) for elem in compound_split if not elem.isdigit()]

                print(f"compound_elements: {compound_elements}")
                #print(f"elems: {elems}")
                if (set(compound_elements) == set(elems)): electrolyte_coeff = coeff
                
                if ((set(compound_elements) == set(['Li'])) 
                    or (set(compound_elements) == set(['Na']))): 
                    print(f"entered open_el_coeff: {coeff}")
                    open_el_coeff = coeff

            #print(f"end of react")
            product_vol = 0
            for compound, coeff in prod_dict.items():
                print(f"product compound: {compound}")
                reduced_compound, gcd = reduce_compound_and_coeff(compound)
                compound_split = seperate_string_number(compound)
                compound_elements = [str(elem) for elem in compound_split if not elem.isdigit()]

                if (reduced_compound not in compound_dict): 
                    KeyError("Key 'compound' is not in dictionary 'compound_dict'")
                    print(f"Key error")
                    exit()
                
                elif ((compound == "Cl2") or (compound == "O2") and (not reservoir_used)): 
                     pass # include this if statement in case of omitting volume of open element
                 
                elif ((set(compound_elements) == set(['Li'])) 
                    or (set(compound_elements) == set(['Na']))): 
                    print(f"entered open_el_coeff: {coeff}")
                    open_el_coeff = -coeff

                else: 
                    molar_volume = float(gcd) * compound_dict[reduced_compound]
                    #print(f"product compound: {compound}")
                    #print(f"product molar_volume: {molar_volume}   product coeff: {coeff}")
                    product_vol += float(molar_volume) * float(coeff)

            print(f"volume: {molar_volume}")                    
            volume_change = product_vol - reactant_vol    
            volume_ratio = product_vol / reactant_vol
            print(f"reactant_vol: {reactant_vol}   product_vol: {product_vol}   volume_change: {volume_change}")
            rel_vol_change_percentage = ( volume_change / reactant_vol ) * 100
            
            react_line_strip = " ".join(react_line.split()[:-2])
            reac_vol_dict[react_line_strip] = rel_vol_change_percentage
            reac_energy_dict[react_line_strip] = react_energy
            reac_electrolyte_mol_dict[react_line_strip] = electrolyte_coeff
            react_open_el_ratio_dict[react_line_strip] = open_el_coeff
            reac_vol_ratio_dict[react_line_strip] = volume_ratio
            vol_by_pot_file.write(react_line_strip + f" volume change: {round(volume_change, 4)}, relative vol change: {round(rel_vol_change_percentage, 4)}% \n ") 

    vol_by_pot_file.write(f"\n\n") 

    return reac_vol_dict, reac_energy_dict, reac_electrolyte_mol_dict, react_open_el_ratio_dict, reac_vol_ratio_dict


def find_max_min_in_dict(dictionary, youngs_mod=[]):
    val_max = 0
    val_min = 0
    ### setting y bounds for all volume plots ###
    for i in range(len(dictionary)):
        for j in range(len(dictionary[0])):
            #print(f"i: {i} j: {j}")
            for reaction, val_change in dictionary[i][j].items():
                if (len(youngs_mod) == 0): pass
                else: val_change = np.cbrt(val_change / 100) * youngs_mod[i]

                if (val_change > val_max): val_max = val_change
                elif (val_change < val_min): val_min = val_change

    return val_min, val_max


def get_x_pos(plt_idx, shape, separation):

    x_pos = plt_idx[0] * shape[1] + ((plt_idx[1] - (shape[1] / 2)) * 0.78) + (shape[1] / 2)

    return x_pos[0]


def create_xlsx_sheet(name, _electrolytes, _voltages):
    # Create a workbook and add a worksheet.
    _workbook = xlsxwriter.Workbook(name)
    _wksheet_energies = _workbook.add_worksheet("reaction energies")
    _wksheet_volumechanges = _workbook.add_worksheet("volume changes")
    _wksheet_strains = _workbook.add_worksheet("strains")
    _wksheet_normstressintensity = _workbook.add_worksheet("strain x modulus")
    _wksheet_volratios = _workbook.add_worksheet("volume ratios")
    _wksheet_reactions = _workbook.add_worksheet("reactions")

    # Add a bold format to use to highlight cells.
    bold = _workbook.add_format({'bold': True})

    # Create a format to use in the merged range.
    merge_format = _workbook.add_format(
        {
            "bold": 1,
            "border": 1,
            "align": "center",
            "valign": "vcenter",
            "rotation": 90
        }
    )

    # Iterate over the data and write it out row by row.
    xlsx_row = 2
    xlsx_col = 1

    _wksheet_energies.merge_range(2, 0, 2+len(_voltages), 0, 'Potential of Li/Li+ or Na/Na+ (eV)', merge_format)
    _wksheet_volumechanges.merge_range(2, 0, 2+len(_voltages), 0, 'Potential of Li/Li+ or Na/Na+ (eV)', merge_format)
    _wksheet_strains.merge_range(2, 0, 2+len(_voltages), 0, 'Potential of Li/Li+ or Na/Na+ (eV)', merge_format)
    _wksheet_normstressintensity.merge_range(2, 0, 2+len(_voltages), 0, 'Potential of Li/Li+ or Na/Na+ (eV)', merge_format)
    _wksheet_volratios.merge_range(2, 0, 2+len(_voltages), 0, 'Potential of Li/Li+ or Na/Na+ (eV)', merge_format)
    _wksheet_reactions.merge_range(2, 0, 2+len(_voltages), 0, 'Potential of Li/Li+ or Na/Na+ (eV)', merge_format)


    # Create a format to use in the merged range.
    merge_format = _workbook.add_format(
        {
            "bold": 1,
            "border": 1,
            "align": "center",
            "valign": "vcenter",
        }
    )


    xlsx_row = 1
    xlsx_col = 2
    _wksheet_energies.merge_range(0, 2, 0, 2+len(_electrolytes),  'Solid electrolytes', merge_format)
    _wksheet_volumechanges.merge_range(0, 2, 0, 2+len(_electrolytes),  'Solid electrolytes', merge_format)
    _wksheet_strains.merge_range(0, 2, 0, 2+len(_electrolytes),  'Solid electrolytes', merge_format)
    _wksheet_normstressintensity.merge_range(0, 2, 0, 2+len(_electrolytes),  'Solid electrolytes', merge_format)
    _wksheet_volratios.merge_range(0, 2, 0, 2+len(_electrolytes),  'Solid electrolytes', merge_format)
    _wksheet_reactions.merge_range(0, 2, 0, 2+len(_electrolytes),  'Solid electrolytes', merge_format)


    xlsx_row = 2
    xlsx_col = 1
    for voltage in (reaction_voltages):
        _wksheet_energies.write(xlsx_row, xlsx_col, voltage, bold)
        _wksheet_volumechanges.write(xlsx_row, xlsx_col, voltage, bold)
        _wksheet_strains.write(xlsx_row, xlsx_col, voltage, bold)
        _wksheet_normstressintensity.write(xlsx_row, xlsx_col, voltage, bold)
        _wksheet_volratios.write(xlsx_row, xlsx_col, voltage, bold)
        _wksheet_reactions.write(xlsx_row, xlsx_col, voltage, bold)
        xlsx_row += 1


    xlsx_row = 1
    xlsx_col = 2
    for electrolyte in _electrolytes:
        # Write some data headers.
        _wksheet_energies.write(xlsx_row, xlsx_col, electrolyte, bold)
        _wksheet_volumechanges.write(xlsx_row, xlsx_col, electrolyte, bold)
        _wksheet_strains.write(xlsx_row, xlsx_col, electrolyte, bold)
        _wksheet_normstressintensity.write(xlsx_row, xlsx_col, electrolyte, bold)
        _wksheet_volratios.write(xlsx_row, xlsx_col, electrolyte, bold)
        _wksheet_reactions.write(xlsx_row, xlsx_col, electrolyte, bold)

    _wksheets = [_wksheet_energies, _wksheet_volumechanges, _wksheet_strains, 
        _wksheet_normstressintensity, _wksheet_volratios, _wksheet_reactions]
    
    return _workbook, _wksheets


### options ###
reservoir_used = True
open_metal_vol_included = True
bare_metals_excluded = True


reaction_elems = np.asarray([set(["Li", "Ge", "P", "S"]), set(["Li", "In", "Cl"]), set(["Li", "La", "Zr", "O"]), 
     set(["Li", "P", "S", "Cl"]), set(["Na", "Sb", "S"]), set(["Na", "Br", "O"]), set(["Al", "Na", "O"]), 
     set(["Li", "P", "S"]), set(["Li", "O", "H", "Cl"]), set(["Na", "Zr", "Si", "O"])])
### set(["Li", "La", "Zr", "Ta", "O"]),
# reaction_elems = np.asarray([set(["Li", "P", "S"])])
### set(["Li", "La", "Zr", "Ta", "O"]),


electrolyte_compounds = np.asarray(["Li10GeP2S12", "Li3InCl6", "Li7La3Zr2O12", 
   "Li6PS5Cl", "Na3SbS4", "Na3BrO", "NaAl11O17", "Li7P3S11", "Li2OHCl", "Na4Zr2(SiO4)3"])
### (1) 21.7 GPa (LGPS),  (2)  44.5, (3) 175.1 (LLZO) ,   (4) 22.1 (LPSCl)  (5) 26.06 (Na3SbS4),  (6)  57.4 (Na3BrO),  (7) 153.94 (NaAl11O17),  (8)  21.9GPa (LPS),   (9) 43.68GPa (Li2OHCl),  


### options for directory for getting reactions (including open metal reservoir) ###
if reservoir_used:

    ### options for directory for plot (including open metal reservoir) ###
    if bare_metals_excluded: 

        ### options for excluding bare metals from reaction products ###
        dirpath = os.getcwd() + "/revised_reactions_nobaremetals_2020compat/"
        bare_metatals_str = "nobaremetals"
        
        if open_metal_vol_included:
            ### including open metal from reactant volume ###
            plot_dirpath = "corrected_selectiveexp_plots_nobaremetals_Li_included_2020compat"
            open_vol_str = "Li_included"
        else:
            ### excluding open metal from reactant volume ###
            plot_dirpath = "corrected_selectiveexp_plots_nobaremetals_Li_excluded_2020compat"
            open_vol_str = "Li_excluded"            
    
    else:

        ### options for including all compounds in reaction products ###
        dirpath = os.getcwd() + "/revised_reactions_theo_exp_2020compat/"
        bare_metatals_str = "theoexp"

        if open_metal_vol_included:
            ### including open metal from reactant volume ###
            plot_dirpath = "corrected_selectiveexp_plots_Li_included_2020compat"
            open_vol_str = "Li_included"
        else:
            ### excluding open metal from reactant volume ###
            plot_dirpath = "corrected_selectiveexp_plots_Li_excluded_2020compat"
            open_vol_str = "Li_excluded"

### options for directory for getting reactions (no open metal reservoir) ###
else:    

    ### options for directory for plot (including open metal reservoir) ###
    if bare_metals_excluded: 
        bare_metatals_str = "nobaremetals"

        ### options for excluding bare metals from reaction products ###
        dirpath = os.getcwd() + "/revised_reactions_nobaremetals_2020compat_noreserve/"
        plot_dirpath = "corrected_selectiveexp_plots_nobaremetals_2020compat_noreserve"        
    
    else:
        bare_metatals_str = "theoexp"

        ### options for including all compounds in reaction products ###
        dirpath = os.getcwd() + "/revised_reactions_theo_exp_2020compat_noreserve/"
        plot_dirpath = "corrected_selectiveexp_plots_2020compat_noreserve"

if reservoir_used:   
    reaction_voltages = np.round(np.arange(-1.0, 1.01, 0.1), 2)
    workbook, wksheets = create_xlsx_sheet(f"electrochemical_embrittlement_{bare_metatals_str}_{open_vol_str}.xlsx", electrolyte_compounds, reaction_voltages)
else:   
        reaction_voltages = np.array([0])
        workbook, wksheets = create_xlsx_sheet(f"electrochemical_embrittlement_{bare_metatals_str}_noReserve.xlsx", electrolyte_compounds, reaction_voltages)

print(plot_dirpath)
# exit()

#plot_dirpath = "cole_LPS_plots"
#dirpath = os.getcwd() + "/cole_LPS_reactions/"

# electrolyte_compounds = np.asarray(["Li7P3S11"])

reaction_voltages = np.round(np.arange(-1.0, 1.01, 0.1), 2)
'''reaction_voltages = np.round(np.logspace(-3,0,10), 3)
reaction_voltages = np.append(0, reaction_voltages)
reaction_voltages = np.concatenate([reaction_voltages, np.round(np.linspace(-4,4,20),3)])
reaction_voltages = np.sort(reaction_voltages)
'''#reaction_voltages = np.asarray([-4, -3.579])
workbook, wksheets = create_xlsx_sheet(f"electrochemical_embrittlement_bare.xlsx", electrolyte_compounds, reaction_voltages)


if os.path.isdir(plot_dirpath): pass
else: os.mkdir(plot_dirpath)

wksheet_energies = wksheets[0]
wksheet_volumechanges = wksheets[1]
wksheet_strains = wksheets[2]
wksheet_normstressintensity = wksheets[3]
wksheet_volratios = wksheets[4]
wksheet_reactions = wksheets[5]

reaction_vol_dicts = np.empty((len(electrolyte_compounds),len(reaction_voltages)),dtype=object)
reaction_energies = np.empty((len(electrolyte_compounds),len(reaction_voltages)),dtype=object)
reaction_elec_mols = np.empty((len(electrolyte_compounds),len(reaction_voltages)),dtype=object)
reaction_open_el_mols = np.empty((len(electrolyte_compounds),len(reaction_voltages)),dtype=object)
reaction_vol_ratios = np.empty((len(electrolyte_compounds),len(reaction_voltages)),dtype=object)

#electrolyte_youngs_mod = np.asarray([37.19, 100, 149.8, 100, 22 , 33.9, 100, 100, 100, 100])
electrolyte_youngs_mod = np.asarray([21.7, 44.5, 175.1, 22.1, 26.06, 57.4, 153.94, 21.9, 43.68, 141.65])


for file in os.listdir(dirpath):

    input_file_split = (file.split("_"))
    input_elems = input_file_split[:-2]
    #sub_dirpath = os.getcwd() + "/exp_volchange_new"
    
    sub_dirpath = os.getcwd() + "/stable_products_volchange_new"
    # sub_dirpath = os.getcwd() + "/test_volume_changes/"

    for sub_file in os.listdir(sub_dirpath):

        print(f"input_elems: {input_elems}")
        print(f"reaction_elems: {reaction_elems}")
        # if (set(input_elems) in reaction_elems):
        sub_file_split = (sub_file.split("_"))
        sub_elems = sub_file_split[:-2]

        if (set(sub_elems) == set(input_elems)):
            
            print(f"sub_elems: {sub_elems}    input_elems: {input_elems}")
            react_vol_dict, react_energy_dict, react_elec_mol_dict, react_open_el_ratio_dict, react_vol_ratio_dict = assign_vol_changes(sub_file, file, dirpath, open_metal_vol_included)
            filename_split = file.split("_")
            elems = set(filename_split[:-2])
            print(f"elems: {elems}")
            
            voltage = float(filename_split[-2].strip("volts"))
            
            if voltage in reaction_voltages:
            
                print(f"voltage: {voltage}")
                elem_i = np.where(reaction_elems == elems)[0][0]
                voltage_i = np.where(reaction_voltages == voltage)[0][0]
                reaction_vol_dicts[elem_i][voltage_i] = react_vol_dict
                reaction_energies[elem_i][voltage_i] = react_energy_dict
                reaction_elec_mols[elem_i][voltage_i] = react_elec_mol_dict
                reaction_open_el_mols[elem_i][voltage_i] = react_open_el_ratio_dict
                reaction_vol_ratios[elem_i][voltage_i] = react_vol_ratio_dict

print(f"reaction_voltages: {reaction_voltages}")
print(f"reaction_vol_dicts: {reaction_vol_dicts}")

### setting y bounds for all volume plots ###
vol_min, vol_max = find_max_min_in_dict(reaction_vol_dicts)

### setting y bounds for stress plots ###
stress_min, stress_max = find_max_min_in_dict(reaction_vol_dicts, youngs_mod=electrolyte_youngs_mod)
energy_min, energy_max = find_max_min_in_dict(reaction_energies)

vol_vals = []
stress_vals = []
energy_vals = []

electrolyte_labels = np.asarray(["Li10GeP2S12", "Li3InCl6", "Li7La3Zr2O12", 
   "Li6PS5Cl", "Na3SbS4", "Na3BrO", "NaAl11O17", "Li7P3S11", "Li2OHCl", "Na4Zr2(SiO4)3"])

grouped_elec_labels = np.asarray( 
    [["Li10GeP2S12", "Li6PS5Cl", "Li7P3S11", "Na3SbS4"], # Sulfides
    ["Li7La3Zr2O12", "Na3BrO", "NaAl11O17", "Na4Zr2(SiO4)3"], # Oxides
    ["Li2OHCl", "Li3InCl6", None, None] # Halides
])


electrolyte_labels = np.asarray(["Li10GeP2S12", "Li3InCl6", "Li7La3Zr2O12", 
   "Li6PS5Cl", "Na3SbS4", "Na3BrO", "NaAl11O17", "Li7P3S11", "Li2OHCl", "Na4Zr2(SiO4)3"])
'''
grouped_elec_labels = np.asarray( 
    [["Li7P3S11"]] # Halides
)

electrolyte_labels = np.asarray(["Li7P3S11"])
'''
subscript_electrolyte_labels = []
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")  
for electrolyte in electrolyte_labels:    
    subscripted_compound = electrolyte.translate(SUB)
    subscript_electrolyte_labels.append(subscripted_compound)


subscript_electrolyte_compounds = []
for electrolyte in electrolyte_compounds:
    subscripted_compound = electrolyte.translate(SUB)
    subscript_electrolyte_compounds.append(subscripted_compound)

# Creating colormap for cumulative plots
cmap_name = "turbo" # "gnuplot2"
cumulative_plot_handles = {}
start_val = 0
stop_val = len(reaction_vol_dicts) + 1
cmap_in = plt.get_cmap(cmap_name)
norm_in = mpl.colors.Normalize(vmin=start_val-1, vmax=stop_val+1)
cumulativeplots_scalarMap = cm.ScalarMappable(norm=norm_in, cmap=cmap_in)

cumul_fig, cumul_ax = plt.subplots((3), figsize=(6,6), constrained_layout = True)
cumul_line_fig, cumul_line_ax = plt.subplots(figsize=(8,6), constrained_layout = True)
trans_strains = np.zeros((len(reaction_vol_dicts), len(reaction_vol_dicts[0])))
vol_ratios = np.zeros((len(reaction_vol_dicts), len(reaction_vol_dicts[0])))

# Creating figure objects for plots for each material
fig, ax = plt.subplots((3), figsize=(10,6), constrained_layout = True)
fontsize_ = 34

# Creating figure objects for plots for each material
reactions_fig, reactions_ax = plt.subplots((2), figsize=(10,6), constrained_layout = True)
#fontsize_ = 56

### voltage for plot of transformation stran, volume change, and reaction energy ###
### spanning all compounds ###
if (reservoir_used): 
    bar_width = reaction_voltages[1] - reaction_voltages[0]
    desired_voltages = [-1, 0, 1]
    x_separation = 0.22
    
else: 
    bar_width = 1
    desired_voltages = [0]
    x_separation = 0.78

hatch_options = ["////", None, "\\\\\\\\"]
_linewidth = 0.5
# desired_voltages = reaction_voltages
# centered_voltage_idx = np.where(reaction_voltages == 0)[0][0]

# hatch_options = ["////"] * len(desired_voltages)
### creating colormaps for groups of electrolyte compounds ###
cmap_name = "Blues" 
start_val = -1
stop_val = len(grouped_elec_labels[0]) + 1
cmap_in = plt.get_cmap(cmap_name)
norm_in = mpl.colors.Normalize(vmin=start_val, vmax=stop_val+1)
blue_scalarMap = cm.ScalarMappable(norm=norm_in, cmap=cmap_in)
blue_rgb_vals = np.array([[123, 174, 232], [57, 136, 227], [95, 104, 227], [39, 51, 227]]) / 256
cmap_name = "Oranges" 
cmap_in = plt.get_cmap(cmap_name)
norm_in = mpl.colors.Normalize(vmin=start_val, vmax=stop_val+1)
orange_scalarMap = cm.ScalarMappable(norm=norm_in, cmap=cmap_in)
green_rgb_vals = np.array([[151, 227, 141], [188, 242, 116], [129, 212, 150], [73, 209, 108]]) / 256
cmap_name = "Greens" 
cmap_in = plt.get_cmap(cmap_name)
norm_in = mpl.colors.Normalize(vmin=start_val, vmax=stop_val+1)
green_scalarMap = cm.ScalarMappable(norm=norm_in, cmap=cmap_in)
orange_rgb_vals = np.array([[235, 95, 52], [235, 125, 52], [235, 155, 52], [235, 185, 52]]) / 256 

custom_colors = False
if not custom_colors: cumul_ax_scalarmaps = [blue_scalarMap, orange_scalarMap, green_scalarMap]
else: cumul_ax_rgbvals = [blue_rgb_vals, orange_rgb_vals, green_rgb_vals]

plot_handles_cumul_allvoltages = []
plot_line_handles_allvoltages = []
cumul_ax_x_ticks = []
cumul_ax_x_ticks_labels = []

for desired_idx in range(len(desired_voltages)):

    plot_handles_cumul = {}
    plot_line_handles = {}

    j_cumulative = np.where(reaction_voltages[::-1] == desired_voltages[desired_idx])[0][0]
    print(f"j_cumulative: {j_cumulative}")
    max_react_string_len = 0
    for i in range(len(reaction_vol_dicts)):
                
        ### creating list (non-repetitive) of all reactions for given solid ###
        ### electrolyte, removing 0-volume change reactions ###
        ax_reactions = []
        ax_energies = np.zeros((len(reaction_vol_dicts[0])))
        ax_electrolyte = " ".join(reaction_elems[i])
        plot_handles = {} #Uncomment for multiple plots
        reaction_plot_handles = {}

        for reaction, volume_change in reaction_vol_dicts[i][j_cumulative].items():
            if (volume_change == 0) or (reaction_elec_mols[i][j_cumulative][str(reaction)] < 1) or (reaction_energies[i][j_cumulative][reaction] > -1e-04):
                pass

            else: 
                ax_reactions.append(str(reaction))
        
        ### removing duplicate reactions with permutations of compounds  ###
        set_reacts = []
        noduplicates_ax_reactions = []
        ax_reactions = np.asarray(list(set(ax_reactions)))
        subscripted_ax_reactions = []

        for j in range(len(ax_reactions)):
            reaction = ax_reactions[j] 
            reactant_dict, product_dict, react_energy = parse_reaction(reaction)
            # is_redox = check_if_redox(reactant_dict, product_dict)
            new_reac_set = set()

            for compound, coeff in product_dict.items():
                new_reac_set.add(compound)
            
            if new_reac_set not in set_reacts: 
                # print(f"ACCEPTED set reaction: {reaction} ")
                set_reacts.append(new_reac_set)
                noduplicates_ax_reactions.append(reaction)
        
        ax_reactions = np.asarray(noduplicates_ax_reactions)

        for react in ax_reactions:        
            subscripted_ax_reactions.append(make_reaction_subscripted(react)) 
        
        print(f"subscripted_ax_reactions: {subscripted_ax_reactions}")
        for react in subscripted_ax_reactions:
            print(f"react: {react}") 
            print(f"len(react): {len(react)}")
            '''
            label = wx.StaticText(panel, label = react, pos = (100,50))    
            dc = wx.ScreenDC()
            size = dc.GetTextExtent(react)
            text_width = size.width
            '''
            if (len(react) > max_react_string_len): 
                max_react_string_len = len(react)
                print(f"max_react_string_len: {max_react_string_len}")


    ### plotting data ###
    for i in range(len(reaction_vol_dicts)):
        plot_idx = np.where((grouped_elec_labels == electrolyte_labels[i]))
        if not custom_colors: cumul_ax_colormap = cumul_ax_scalarmaps[plot_idx[0][0]]
        else: cumul_ax_colormap = cumul_ax_rgbvals[plot_idx[0][0]]
        x_pos = get_x_pos(plot_idx, grouped_elec_labels.shape, x_separation)

        if (electrolyte_labels[i].translate(SUB) not in cumul_ax_x_ticks_labels):
            subscripted_compound = electrolyte_labels[i].translate(SUB)
            cumul_ax_x_ticks.append(x_pos)
            cumul_ax_x_ticks_labels.append(subscripted_compound)

        print(f"plot_idx: {plot_idx}")
        print(reaction_vol_dicts.size)
        print(f"\n\ni: {i}  electrolyte_compounds[i]: {electrolyte_compounds[i]}")
        
        ### creating list (non-repetitive) of all reactions for given solid ###
        ### electrolyte, removing 0-volume change reactions ###
        ax_reactions = []
        ax_energies = np.zeros((len(reaction_vol_dicts[0])))
        ax_open_el = np.zeros((len(reaction_vol_dicts[0])))
        ax_electrolyte_mols = np.zeros((len(reaction_vol_dicts[0])))
        ax_electrolyte = " ".join(reaction_elems[i])
        plot_handles = {} #Uncomment for multiple plots

        for j in range(len(reaction_vol_dicts[0])):
            print(f"first voltages: {reaction_voltages[j]}")
            for reaction, volume_change in reaction_vol_dicts[i][j].items():

                if (volume_change == 0) or (reaction_elec_mols[i][j][str(reaction)] < 1) or (reaction_energies[i][j][reaction] > -1e-04):
                    # print(f"pass")
                    pass

                else: 
                    # print(f"ACCEPTED reaction: {reaction}  volume_change: {volume_change}   reaction_energies[i][j][str(reaction)]: {reaction_energies[i][j][str(reaction)]}")
                    ax_reactions.append(str(reaction))
                    ax_energies[j] = (reaction_energies[i][j][str(reaction)])
                    ax_open_el[j] = (reaction_open_el_mols[i][j][str(reaction)])
                    ax_electrolyte_mols[j] = (reaction_elec_mols[i][j][str(reaction)]) 
        
        print(f"ax_energies: {ax_energies}")
        ### removing duplicate reactions with permutations of compounds  ###
        set_reacts = []
        noduplicates_ax_reactions = []
        ax_reactions = np.asarray(list(set(ax_reactions)))
        subscripted_ax_reactions = []

        for j in range(len(ax_reactions)):
            reaction = ax_reactions[j] 
            reactant_dict, product_dict, react_energy = parse_reaction(reaction)
            # is_redox = check_if_redox(reactant_dict, product_dict)
            new_reac_set = set()

            for compound, coeff in product_dict.items():
                new_reac_set.add(compound)
            
            if new_reac_set not in set_reacts: 
                # print(f"ACCEPTED set reaction: {reaction} ")
                set_reacts.append(new_reac_set)
                noduplicates_ax_reactions.append(reaction)
        
        ax_reactions = np.asarray(noduplicates_ax_reactions)

        for react in ax_reactions:        
            subscripted_ax_reactions.append(make_reaction_subscripted(react)) 
        
        '''
        for react in subscripted_ax_reactions:
            if (len(react) > max_react_string_len): 
                print(f"react: {react}")
                print(f"len(react): {len(react)}")
                max_react_string_len = len(react)
                print(f"max_react_string_len: {max_react_string_len}")
        '''
        
        # uncomment for multiple plots
        start_val = -1
        stop_val = len(subscripted_ax_reactions)
        cmap_in = plt.get_cmap("turbo")
        norm_in = mpl.colors.Normalize(vmin=start_val-1, vmax=stop_val+1)
        scalarMap = cm.ScalarMappable(norm=norm_in, cmap=cmap_in)

        ### getting minimum energy reactions for cumulative plot ###
        vol_changes_unsorted = []
        reacts_unsorted = []
        energies_unsorted = []
        elec_mols_unsorted = []
        open_el_mols_unsorted = []
        vol_ratios_unsorted = []
        
        ### removing 0-volume reactions and selecting min-energy reaction ###
        ### at each given voltage (multiple can occur) ###
        if (len(reaction_vol_dicts[i][j_cumulative]) > 1):
            min_energy = sys.maxsize
            min_e_reaction = ""
            min_e_vol_change = 0
            min_e_elec_mols = 0
            
            for reaction, volume_change in reaction_vol_dicts[i][j_cumulative].items():
                print(f"reaction_energies[i][j_cumulative][reaction]: {reaction_energies[i][j_cumulative][reaction]}")
                print(f"react_vol_ratio_dict: {react_vol_ratio_dict}")
                if (volume_change == 0) or (reaction_elec_mols[i][j_cumulative][reaction] < 1) or (reaction_energies[i][j_cumulative][reaction] > -1e-04): 
                    pass
                elif (reaction_energies[i][j_cumulative][reaction] < min_energy):
                    min_energy = reaction_energies[i][j_cumulative][reaction]
                    min_e_reaction = reaction
                    min_e_vol_change = volume_change
                    min_e_elec_mols = reaction_elec_mols[i][j_cumulative][reaction]
                    min_e_vol_ratio = reaction_vol_ratios[i][j_cumulative][reaction]

            if not min_e_reaction:
                print(f"in not min e reaction")
                vol_changes_unsorted.append(0)
                reacts_unsorted.append("None")      
                energies_unsorted.append(0)
                elec_mols_unsorted.append(1)
                vol_ratios_unsorted.append(0)
            else:
                vol_changes_unsorted.append(min_e_vol_change)
                reacts_unsorted.append(min_e_reaction)      
                energies_unsorted.append(min_energy)
                elec_mols_unsorted.append(min_e_elec_mols)
                vol_ratios_unsorted.append(min_e_vol_ratio)
        
        else:
            for reaction, volume_change in reaction_vol_dicts[i][j_cumulative].items():
                
                if (volume_change == 0) or (reaction_elec_mols[i][j_cumulative][reaction] != 1) or (reaction_energies[i][j_cumulative][reaction] > -1e-04):
                    vol_changes_unsorted.append(0)
                    reacts_unsorted.append("None")
                    energies_unsorted.append(0)
                    elec_mols_unsorted.append(0)
                    vol_ratios_unsorted.append(0)
                else: 
                    vol_changes_unsorted.append(volume_change)
                    reacts_unsorted.append(reaction)
                    energies_unsorted.append(reaction_energies[i][j_cumulative][reaction])
                    elec_mols_unsorted.append(reaction_elec_mols[i][j_cumulative][reaction])
                    print(f"reaction: {reaction}")
                    vol_ratios_unsorted.append(reaction_vol_ratios[i][j_cumulative][reaction])
        
        reacts_sorted_cumul = [react for _, react in sorted(zip(vol_changes_unsorted, reacts_unsorted))]
        print(f"reacts_sorted_cumul: {reacts_sorted_cumul}")
        energies_sorted_cumul = [react for _, react in sorted(zip(vol_changes_unsorted, energies_unsorted))]
        elec_mols_sorted_cumul = [react for _, react in sorted(zip(vol_changes_unsorted, elec_mols_unsorted))]
        scaled_volumes_sorted_cumul = np.sort(vol_changes_unsorted) 

        ### creating bars ###
        for j in range(len(reaction_vol_dicts[0])):  ### TAB BACK THE REST OF THIS FOR LOOP
            j_reversed = copy.deepcopy(int(len(reaction_voltages) - 1 - j))
            ax_voltage = reaction_voltages[int(len(reaction_voltages) - 1 - j)] 
            print(f"reaction_vol_dicts[0]: {reaction_vol_dicts[0]}")
            print(f"ax_voltage: {ax_voltage}")
            print(f"second voltages: {reaction_voltages[j_reversed]}")
            print(f"cumulative voltage: {reaction_voltages[j_cumulative]}")

            color_i = 0
            
            vol_changes_unsorted = []
            reacts_unsorted = []
            energies_unsorted = []
            elec_mols_unsorted = []
            open_el_mols_unsorted = []
            vol_ratios_unsorted = []
            
            ### removing 0-volume reactions and selecting min-energy reaction ###
            ### at each given voltage (multiple can occur) ###
            if (len(reaction_vol_dicts[i][j_reversed]) > 1):
                min_energy = sys.maxsize
                min_e_reaction = ""
                min_e_vol_change = 0
                min_e_elec_mols = 0
                min_e_open_el_mols = 0
                print(f"statement 1")
                
                for reaction, volume_change in reaction_vol_dicts[i][j_reversed].items():
                    # print(f"reaction: {reaction}  reaction_energies[i][j][reaction]: {reaction_energies[i][j][reaction]}")
                    if (volume_change == 0) or (reaction_elec_mols[i][j_reversed][reaction] < 1) or (reaction_energies[i][j_reversed][reaction] > -1e-04): 
                        print(f"remove energy: {reaction_energies[i][j_reversed][reaction]}")
                        print(f"remove reaction: {reaction}")
                        pass
                    elif (reaction_energies[i][j_reversed][reaction] < min_energy):
                        min_energy = reaction_energies[i][j_reversed][reaction]
                        min_e_reaction = reaction
                        min_e_vol_change = volume_change
                        min_e_elec_mols = reaction_elec_mols[i][j_reversed][reaction]
                        min_e_open_el_mols = reaction_open_el_mols[i][j_reversed][reaction]
                        min_e_vol_ratio = reaction_vol_ratios[i][j_reversed][reaction]
                
                print(f"min_e_reaction: {min_e_reaction}")
                print(f"min_energy: {min_energy}")
                print(f"min_e_vol_change: {min_e_vol_change}")
                print(f"min_e_elec_mols: {min_e_elec_mols}")
                
                if not min_e_reaction: 
                    print(f"in not min e reaction")
                    vol_changes_unsorted.append(0)
                    reacts_unsorted.append("None")      
                    energies_unsorted.append(0)
                    elec_mols_unsorted.append(1)
                    vol_ratios_unsorted.append(0)
                    
                else:
                    print(f"found min e reaction")
                    vol_changes_unsorted.append(min_e_vol_change)
                    reacts_unsorted.append(min_e_reaction)      
                    energies_unsorted.append(min_energy)
                    elec_mols_unsorted.append(min_e_elec_mols)
                    vol_ratios_unsorted.append(min_e_vol_ratio)
            
            else:
                print(f"statement 2")
                print(f"min_e_reaction: {min_e_reaction}")
                print(f"min_energy: {min_energy}")
                print(f"min_e_vol_change: {min_e_vol_change}")
                print(f"min_e_elec_mols: {min_e_elec_mols}")
                for reaction, volume_change in reaction_vol_dicts[i][j_reversed].items():
                    if (volume_change == 0) or (reaction_elec_mols[i][j_reversed][reaction] != 1) or (reaction_energies[i][j_reversed][reaction] > -1e-04):
                        vol_changes_unsorted.append(0)
                        reacts_unsorted.append("None")
                        energies_unsorted.append(0)
                        elec_mols_unsorted.append(0)
                        open_el_mols_unsorted.append(0)
                        vol_ratios_unsorted.append(0)
                        
                    else: 
                        vol_changes_unsorted.append(volume_change)
                        reacts_unsorted.append(reaction)
                        energies_unsorted.append(reaction_energies[i][j_reversed][reaction])
                        elec_mols_unsorted.append(reaction_elec_mols[i][j_reversed][reaction])
                        open_el_mols_unsorted.append(reaction_open_el_mols[i][j_reversed][reaction])
                        print(f"reaction: {reaction}")
                        min_e_vol_ratio = reaction_vol_ratios[i][j_reversed][reaction]
                        vol_ratios_unsorted.append(min_e_vol_ratio)
            
            
            reacts_sorted = [react for _, react in sorted(zip(vol_changes_unsorted, reacts_unsorted))]
            print(f"reacts_sorted: {reacts_sorted}")
            energies_sorted = [react for _, react in sorted(zip(vol_changes_unsorted, energies_unsorted))]
            elec_mols_sorted = [react for _, react in sorted(zip(vol_changes_unsorted, elec_mols_unsorted))]
            open_el_mols_sorted = [react for _, react in sorted(zip(vol_changes_unsorted, open_el_mols_unsorted))]
            scaled_vol_ratios_sorted = [react for _, react in sorted(zip(vol_changes_unsorted, vol_ratios_unsorted))]
            scaled_volumes_sorted = np.sort(vol_changes_unsorted) 

            ### (deprecated) if more than one reaction at given voltage, scaling ###
            ###   volume changes according to abundance based on energies of reactions ###
            T = 300
            kB = 8.6173e-5
            total_vol_change = 0
            total_stress = 0


            ### plotting bars according to volume change for each reaction ###
            if (not reacts_sorted):            
                space_shift = max_react_string_len - len(subscripted_ax_reactions[color_i]) + 3
                labelval = f"{subscripted_ax_reactions[color_i]};"
                for space in range(3): labelval += " "
                labelval += f"{np.round(ax_energies[j_cumulative],3)} eV"
                print(f"len(labelval): {len(labelval)}")
                bar_handle = ax[0].bar(x_pos, 0, width = 1, color = scalarMap.to_rgba(j), alpha = 1, label=f"{labelval}")
                if not custom_colors: bar_handle = cumul_ax[0].bar(x_pos + desired_voltages[desired_idx]*x_separation, 0, width = x_separation, color = cumul_ax_colormap.to_rgba(plot_idx[1][0]), hatch = hatch_options[desired_idx], alpha = 1, label=f"{labelval}", edgecolor = "black", linewidth = _linewidth)
                else: bar_handle = cumul_ax[0].bar(x_pos + desired_voltages[desired_idx]*x_separation, 0, width = x_separation, color = cumul_ax_colormap[plot_idx[1][0]], hatch = hatch_options[desired_idx], alpha = 1, label=f"{labelval}", edgecolor = "black", linewidth = _linewidth)
                vert_shift = 0
                wksheet_energies.write((j+2), (i+2), "None")
                print(f"react sorted empty")
                
            for k in range(len(reacts_sorted)):            
                if (reacts_sorted[k] == "None"):
                    print(f"statement 1 reacts_sorted[k] == None ")
                    react_energy = 0
                    product_dict = {}
                    reactant_dict = {}
                else: 
                    reactant_dict, product_dict, react_energy = parse_reaction(reacts_sorted[k])
                
                new_reac_set = set()
                
                for compound, coeff in product_dict.items():
                    new_reac_set.add(compound)

                if (reacts_sorted[k] == "None"):
                    print(f"reacts_sorted[k] == None statement 1")
                    color_i = 0
                else: 
                    color_i = np.where(np.asarray(set_reacts) == new_reac_set)[0][0]
                
                if (reacts_sorted[k] == "None"):
                    print(f"None statement 3 bar {i}")
                    
                    labelval = f"No reaction"
                    
                    deltaK_div_thickness =  0
                    bar_handle = ax[2].bar(ax_voltage, deltaK_div_thickness, width = round(bar_width, 2), color = scalarMap.to_rgba(color_i), alpha = 1, label=labelval)
                    wksheet_normstressintensity.write((j+2), (i+2), deltaK_div_thickness)

                    labelval = f"No reaction"
                    bar_handle = ax[1].bar(ax_voltage, scaled_volumes_sorted[k], width = round(bar_width, 2), color = scalarMap.to_rgba(color_i), alpha = 1, label=f"{labelval}")
                                        
                    bar_handle = ax[0].bar(ax_voltage, 0, width = round(bar_width, 2), color = scalarMap.to_rgba(color_i), alpha = 1, label=f"{labelval}")
                    bar_line_handle_cumul = cumul_line_ax.scatter(reaction_voltages[j], 0, color = cumulativeplots_scalarMap.to_rgba(i), alpha = 1, label=f"{labelval}", s=8, marker='s')
                    print(f"scaled_volumes_sorted: {scaled_volumes_sorted}")
                    trans_strains[i][j] = 0
                    vol_ratios[i][j] = 1

                    wksheet_energies.write((j+2), (i+2), "None")
                    wksheet_volumechanges.write((j+2), (i+2), "None")
                    wksheet_strains.write((j+2), (i+2), "None")
                    wksheet_volratios.write((j+2), (i+2), "None")
                    wksheet_reactions.write((j+2), (i+2), "None")

                elif (scaled_volumes_sorted[k] == 0): 
                    print(f"error: vol_changes_sorted[k] should have no 0 elements")
                    pass
                else:
                    print(f"reactions not empty")
                    
                    space_shift = max_react_string_len - len(subscripted_ax_reactions[color_i]) + 3
                    labelval = f"{subscripted_ax_reactions[color_i]};"
                    for space in range(3): labelval += " "
                    labelval += f"{np.round(ax_energies[j_cumulative],3)} eV"

                    bar_handle = ax[1].bar(ax_voltage, scaled_volumes_sorted[k], width = round(bar_width, 2), color = scalarMap.to_rgba(color_i), alpha = 1, label=f"{labelval}")
                    
                    if (bar_handle not in plot_handles): plot_handles[labelval] = bar_handle
                    print(f"bar_handle: {bar_handle}")
                    

                    num_open_el_mols = ax_open_el[j_reversed]
                    num_elec_mols = ax_electrolyte_mols[j_reversed]
                    elec_formula_split = seperate_string_number(electrolyte_compounds[i])
                    print(f"elec_formula_split: {elec_formula_split}")

                    elec_el_count = 0
                    el_idx = 0
                    reaction_energies_shifted = []

                    while (el_idx < len(elec_formula_split)):
                        elem = elec_formula_split[el_idx]
                        if elem in set(['Li', 'Na']): el_idx += 1 
                        elif elem.isdigit(): elec_el_count += int(elem)
                        el_idx += 1 
                    print(f"j: {j}    ax_energies[j]: {ax_energies[j]} ") 
                    for mu in reaction_voltages:
                        reaction_energy_shifted = ax_energies[j_reversed] - (num_open_el_mols * (mu-ax_voltage)) / ( num_elec_mols * elec_el_count)
                        reaction_energies_shifted.append(reaction_energy_shifted)
                        print(f"mu: {mu}    reaction_energy_shifted: {reaction_energy_shifted}  num_open_el_mols: {num_open_el_mols}")
                    
                    print(f"ax_voltage: {ax_voltage}")
                    print(f"scaled_volumes_sorted[k]: {scaled_volumes_sorted[k]}")
                    shifted_reaction_labelval = f"{subscripted_ax_reactions[color_i]}"
                    reactions_bar_handle = reactions_ax[0].plot(reaction_voltages, reaction_energies_shifted, color = scalarMap.to_rgba(color_i), label=f"{shifted_reaction_labelval}")                    
                    if (reactions_bar_handle[0] not in reaction_plot_handles): reaction_plot_handles[shifted_reaction_labelval] = reactions_bar_handle[0]
                    
                    # if (j == 0): bar_width = reaction_voltages[1] - reaction_voltages[0] 
                    # elif (j == (len(reaction_voltages)-1)): bar_width = reaction_voltages[-1] - reaction_voltages[-2] 
                    # else: bar_width = reaction_voltages[j] - reaction_voltages[j+1] 
                    
                    reactions_bar_handle = reactions_ax[1].scatter(ax_voltage, scaled_volumes_sorted[k], color = scalarMap.to_rgba(color_i), label=f"{shifted_reaction_labelval}") 
                    print(f"reactions_bar_handle: {reactions_bar_handle}")
                    # if (bar_handle not in plot_handles): plot_handles[labelval] = bar_handle
                    
                    total_vol_change += scaled_volumes_sorted[k]
                    bar_handle = ax[0].bar(ax_voltage, np.round(ax_energies[j_reversed],3), width = round(bar_width, 2), color = scalarMap.to_rgba(color_i), alpha = 1, label=f"{labelval}")
                    bar_line_handle_cumul = cumul_line_ax.scatter(reaction_voltages[j], np.cbrt(scaled_volumes_sorted[k] / 100), color = cumulativeplots_scalarMap.to_rgba(i), alpha = 1, label=f"{labelval}", s=8, marker='s')
                    print(f"scaled_volumes_sorted: {scaled_volumes_sorted}")
                    trans_strains[i][j] = np.cbrt(scaled_volumes_sorted[k] / 100)
                    vol_ratios[i][j] = scaled_vol_ratios_sorted[k]

                    wksheet_energies.write((j+2), (i+2), np.round(ax_energies[j_reversed],3))
                    wksheet_volumechanges.write((j+2), (i+2), scaled_volumes_sorted[k])
                    wksheet_strains.write((j+2), (i+2), trans_strains[i][j_reversed])
                    wksheet_volratios.write((j+2), (i+2), vol_ratios[i][j_reversed])
                    wksheet_reactions.write((j+2), (i+2), labelval)

                    if (bar_line_handle_cumul not in plot_line_handles): plot_line_handles[subscript_electrolyte_compounds[i]] = bar_line_handle_cumul

            for k in range(len(reacts_sorted_cumul)):

                print(f"k: {k}, reacts_sorted_cumul: {reacts_sorted_cumul}")

                if (reacts_sorted_cumul[k] == "None"):
                    color_i = 0
                    react_energy = 0
                    product_dict = {}
                    reactant_dict = {}
                else:
                    reactant_dict, product_dict, react_energy = parse_reaction(reacts_sorted_cumul[k])
                    new_reac_set = set()
                    for compound, coeff in product_dict.items():
                        new_reac_set.add(compound)

                    color_i = np.where(np.asarray(set_reacts) == new_reac_set)[0][0]
                            
                if (reacts_sorted_cumul[k] == "None"):
                    print(f"None statement 2 bar for i: {i}")
                    print(f"None statement 2 x_pos: {x_pos}")
                    print(f"None statement 2 scaled_volumes_sorted_cumul[k]: {scaled_volumes_sorted_cumul[k]}")
                    labelval = f"No {subscript_electrolyte_compounds[i]} reaction"
                    if not custom_colors: bar_handle_cumul = cumul_ax[1].bar(x_pos + desired_voltages[desired_idx]*x_separation, 0, width = x_separation, color = cumul_ax_colormap.to_rgba(plot_idx[1][0]), hatch = hatch_options[desired_idx], alpha = 1, label=labelval, edgecolor = "black", linewidth = _linewidth)
                    else: bar_handle_cumul = cumul_ax[1].bar(x_pos + desired_voltages[desired_idx]*x_separation, 0, width = x_separation, color = cumul_ax_colormap[plot_idx[1][0]], hatch = hatch_options[desired_idx], alpha = 1, label=labelval, edgecolor = "black", linewidth = _linewidth)
                    

                    if (bar_handle_cumul not in plot_handles_cumul): plot_handles_cumul[labelval] = bar_handle_cumul

                    total_vol_change += scaled_volumes_sorted_cumul[k]
                    if not custom_colors: bar_handle_cumul = cumul_ax[0].bar(x_pos + desired_voltages[desired_idx]*x_separation, 0, width = x_separation, color = cumul_ax_colormap.to_rgba(plot_idx[1][0]), hatch = hatch_options[desired_idx], alpha = 1, label=labelval, edgecolor = "black", linewidth = _linewidth)
                    else: bar_handle_cumul = cumul_ax[0].bar(x_pos + desired_voltages[desired_idx]*x_separation, 0, width = x_separation, color = cumul_ax_colormap[plot_idx[1][0]], hatch = hatch_options[desired_idx], alpha = 1, label=labelval, edgecolor = "black", linewidth = _linewidth)
                
                elif (scaled_volumes_sorted_cumul[k] == 0): 
                    print(f"error: vol_changes_sorted[k] should have no 0 elements")
                    pass
                else:
                    print(f"bar for i: {i}")
                    print(f"x_pos: {x_pos}")
                    print(f"scaled_volumes_sorted_cumul[k]: {scaled_volumes_sorted_cumul[k]}")   
                    print(f"ax_energies[j_cumulative]: {ax_energies[j_cumulative]}")    
                    space_shift = max_react_string_len - len(subscripted_ax_reactions[color_i]) + 3
                    labelval = f"{subscripted_ax_reactions[color_i]};"
                    for space in range(3): labelval += " "
                    labelval += f"{np.round(ax_energies[j_cumulative],3)} eV" 
                    print(f"labelval: {labelval}")    
                    print(f"len(labelval): {len(labelval)}")  
                    if not custom_colors: bar_handle_cumul = cumul_ax[1].bar(x_pos + desired_voltages[desired_idx]*x_separation, scaled_volumes_sorted_cumul[k], width = x_separation, color = cumul_ax_colormap.to_rgba(plot_idx[1][0]), hatch = hatch_options[desired_idx], alpha = 1, label=f"{labelval}", edgecolor = "black", linewidth = _linewidth)
                    else: bar_handle_cumul = cumul_ax[1].bar(x_pos + desired_voltages[desired_idx]*x_separation, scaled_volumes_sorted_cumul[k], width = x_separation, color = cumul_ax_colormap[plot_idx[1][0]], hatch = hatch_options[desired_idx], alpha = 1, label=f"{labelval}", edgecolor = "black", linewidth = _linewidth)
                    

                    if (scaled_volumes_sorted_cumul > 0): vert_shift = scaled_volumes_sorted_cumul[k] + 2.5
                    else: vert_shift = scaled_volumes_sorted_cumul[k] - 12.5
                    # cumul_ax[1].text(x_pos + desired_voltages[desired_idx]*x_separation, vert_shift, f"{int(desired_voltages[desired_idx])}V", ha="center", va="bottom", rotation=60, fontsize= 6) #, bbox=dict(facecolor='white', alpha=0.5))
                    vol_vals.append(scaled_volumes_sorted_cumul[k])

                    print(f"x_pos + desired_voltages[desired_idx]*x_separation: {x_pos + desired_voltages[desired_idx]*x_separation}")
                    if (bar_handle_cumul not in plot_handles_cumul): plot_handles_cumul[labelval] = bar_handle_cumul

                    total_vol_change += scaled_volumes_sorted_cumul[k]
                    if not custom_colors: bar_handle_cumul = cumul_ax[0].bar(x_pos + desired_voltages[desired_idx]*x_separation, np.round(ax_energies[j_cumulative],3), width = x_separation, color = cumul_ax_colormap.to_rgba(plot_idx[1][0]), hatch = hatch_options[desired_idx], alpha = 1, label=f"{labelval}", edgecolor = "black", linewidth = _linewidth)
                    else: bar_handle_cumul = cumul_ax[0].bar(x_pos + desired_voltages[desired_idx]*x_separation, np.round(ax_energies[j_cumulative],3), width = x_separation, color = cumul_ax_colormap[plot_idx[1][0]], hatch = hatch_options[desired_idx], alpha = 1, label=f"{labelval}", edgecolor = "black", linewidth = _linewidth)
                    
                    if (ax_energies[j_cumulative] > 0): vert_shift = ax_energies[j_cumulative] + 0.6
                    else: vert_shift = vert_shift = ax_energies[j_cumulative] - 0.6
                    # cumul_ax[0].text(x_pos + desired_voltages[desired_idx]*x_separation, vert_shift, f"{int(desired_voltages[desired_idx])}V", ha="center", va="bottom", rotation=60, fontsize= 6) #, bbox=dict(facecolor='white', alpha=0.5))
                    energy_vals.append(scaled_volumes_sorted_cumul[k])

            print(f"plot handles cumulative mid: {plot_handles_cumul}")

            for k in range(len(reacts_sorted)):

                print(f"statement 3 k: {k}, reacts_sorted: {reacts_sorted}")
                if (reacts_sorted[k] == "None"): 
                    print(f"None statement 3 idx {i}")
                    color_i = 0
                    react_energy = 0
                    product_dict = {}
                    reactant_dict = {}                
                else:
                    reactant_dict, product_dict, react_energy = parse_reaction(reacts_sorted[k])
                
                new_reac_set = set()

                for compound, coeff in product_dict.items():
                    new_reac_set.add(compound)

                if (reacts_sorted[k] == "None"): color_i = 0
                else: color_i = np.where(np.asarray(set_reacts) == new_reac_set)[0][0]

                if (reacts_sorted[k] == "None"):
                    print(f"None statement 3 bar {i}")
                    
                    labelval = f"No reaction"
                    
                    deltaK_div_thickness =  0
                    bar_handle = ax[2].bar(ax_voltage, deltaK_div_thickness, width = round(bar_width, 2), color = scalarMap.to_rgba(color_i), alpha = 1, label=labelval)
                    wksheet_normstressintensity.write((j+2), (i+2), "None")

                elif (reacts_sorted[k] == 0): 
                    print(f"error: vol_changes_sorted[k] should have no 0 elements")
                    pass
                else:
                    ###  calculating change in effective fracture toughness  ###
                    ### redox active ==> reaction at tip, chemical ==> reaction along crack length ###
                    space_shift = max_react_string_len - len(subscripted_ax_reactions[color_i]) + 3
                    labelval = f"{subscripted_ax_reactions[color_i]};"
                    for space in range(3): labelval += " "
                    labelval += f"{np.round(ax_energies[j_cumulative],3)} eV"
                    #alpha = 0.42
                    #thickness = 3.75e-8
                    deltaK_div_thickness =  np.cbrt(scaled_volumes_sorted[k] / 100) * electrolyte_youngs_mod[i]
                    bar_handle = ax[2].bar(ax_voltage, deltaK_div_thickness, width = round(bar_width, 2), color = scalarMap.to_rgba(color_i), alpha = 1, label=f"{ax_reactions[color_i]}")
                    wksheet_normstressintensity.write((j+2), (i+2), deltaK_div_thickness)

            for k in range(len(reacts_sorted_cumul)):
                print(f"statement 4 bar for i: {i}")
                print(f"statement 4 x_pos: {x_pos}")
                print(f"statement 4 scaled_volumes_sorted_cumul[k]: {scaled_volumes_sorted_cumul[k]}")

                print(f"k: {k}, reacts_sorted_cumul: {reacts_sorted_cumul}")

                if (reacts_sorted_cumul[k] == "None"): 
                    print(f"None statement 3 idx {i}")
                    color_i = 0
                    react_energy = 0
                    product_dict = {}
                    reactant_dict = {}  

                else:
                    reactant_dict, product_dict, react_energy = parse_reaction(reacts_sorted_cumul[k])
                
                new_reac_set = set()
                for compound, coeff in product_dict.items():
                    new_reac_set.add(compound)

                
                if (reacts_sorted_cumul[k] == "None"):
                    color_i = 0

                else:
                    color_i = np.where(np.asarray(set_reacts) == new_reac_set)[0][0]

                if (reacts_sorted_cumul[k] == "None"):
                    print(f"None statement 4 bar {i}")
                    
                    labelval = f"No {subscript_electrolyte_compounds[i]} reaction"
                    deltaK_div_thickness = 0
                    if not custom_colors: bar_handle_cumul = cumul_ax[2].bar(x_pos + desired_voltages[desired_idx]*x_separation, deltaK_div_thickness, width = x_separation, color = cumul_ax_colormap.to_rgba(plot_idx[1][0]), hatch = hatch_options[desired_idx], alpha = 1, label=labelval, edgecolor = "black", linewidth = _linewidth)           
                    else: bar_handle_cumul = cumul_ax[2].bar(x_pos + desired_voltages[desired_idx]*x_separation, deltaK_div_thickness, width = x_separation, color = cumul_ax_colormap[plot_idx[1][0]], hatch = hatch_options[desired_idx], alpha = 1, label=labelval, edgecolor = "black", linewidth = _linewidth)
                    
                elif (scaled_volumes_sorted_cumul[k] == 0): 
                    print(f"error: vol_changes_sorted[k] should have no 0 elements")
                    pass
                
                else:
                    ###  calculating change in effective fracture toughness  ###
                    ### redox active ==> reaction at tip, chemical ==> reaction along crack length ###
                    space_shift = max_react_string_len - len(subscripted_ax_reactions[color_i]) + 3
                    labelval = f"{subscripted_ax_reactions[color_i]};"
                    for space in range(3): labelval += " "
                    labelval += f"{np.round(ax_energies[j_cumulative],3)} eV"
                    print(f"len(labelval): {len(labelval)}")
                    #alpha = 0.42
                    #thickness = 3.75e-8
                    deltaK_div_thickness =  np.cbrt(scaled_volumes_sorted_cumul[k] / 100) * electrolyte_youngs_mod[i]
                    if not custom_colors: bar_handle_cumul = cumul_ax[2].bar(x_pos + desired_voltages[desired_idx]*x_separation, deltaK_div_thickness, width = x_separation, color = cumul_ax_colormap.to_rgba(plot_idx[1][0]), hatch = hatch_options[desired_idx], alpha = 1, label=f"{ax_reactions[color_i]}", edgecolor = "black", linewidth = _linewidth)
                    else: bar_handle_cumul = cumul_ax[2].bar(x_pos + desired_voltages[desired_idx]*x_separation, deltaK_div_thickness, width = x_separation, color = cumul_ax_colormap[plot_idx[1][0]], hatch = hatch_options[desired_idx], alpha = 1, label=f"{ax_reactions[color_i]}", edgecolor = "black", linewidth = _linewidth)
                    
                    if (deltaK_div_thickness > 0): vert_shift = deltaK_div_thickness + 10
                    else: vert_shift = deltaK_div_thickness - 22.5
                    # cumul_ax[2].text(x_pos + desired_voltages[desired_idx]*x_separation, vert_shift, f"{int(desired_voltages[desired_idx])}V", ha="center", va="bottom", rotation=60, fontsize= 6) #, bbox=dict(facecolor='white', alpha=0.25,))
                    stress_vals.append(deltaK_div_thickness)

            ### setting yscale for axis object ###
            if (total_vol_change > vol_max): vol_max = total_vol_change
            elif (total_vol_change < vol_min): vol_min = total_vol_change

            if (total_stress > stress_max): stress_max = total_stress   
            elif (total_stress < stress_min): stress_min = total_stress
            print(f"plot handles cumulative: {plot_handles_cumul}")

        print("\n\n")
        fontsize_ = 12
        plt.rc('font', size=fontsize_) 

        ### Plotting reactions for indivdual materials ###
        ax[0].set_title(subscript_electrolyte_labels[i])
        #ax[0].set_xlim(np.min(reaction_voltages)-0.25,np.max(reaction_voltages)+0.25)
        ax[0].set_ylabel("Reaction \nEnergy \n(eV/atom)", size = fontsize_)

        #ax[1].set_xlim(np.min(reaction_voltages)-0.25,np.max(reaction_voltages)+0.25)
        ax[1].set_ylabel("Relative \nvolume \nchange (%)", size = fontsize_)
        
        #ax[2].set_xlim(np.min(reaction_voltages)-0.25,np.max(reaction_voltages)+0.25)
        ax[2].set_ylabel(r"$ε_{trans}E_{youngs}$ (GPa)", size = fontsize_)
        ax[2].set_xlabel("Voltage of Li or Na relative to Li/Li+ or Na/Na+ (V)")
        
        figure = plt.gcf() # get current figure
        figure.set_size_inches(10, 6)

        fig.savefig(f"{plot_dirpath}/stress_electrolyte_{electrolyte_compounds[i]}.png",dpi=300)

        legend = ax[1].legend(plot_handles.values(), plot_handles.keys(), bbox_to_anchor=(-0.15, 0), loc="lower left", fontsize = "xx-small")
        legend.get_frame().set_alpha(1)            
        figure = plt.gcf() # get current figure
        fig.savefig(f"{plot_dirpath}/legend_volchange_electrolyte_{electrolyte_compounds[i]}.png",dpi=300)
        
        ax[0].clear()
        ax[1].clear()
        ax[2].clear()
    
    plot_handles_cumul_allvoltages.append(plot_handles_cumul)
    plot_line_handles_allvoltages.append(plot_line_handles)

### Plotting cumulative plots ###
electrolyte_labels = np.asarray(["Li10GeP2S12", "Li3InCl6", "Li7La3Zr2O12", 
   "Li6PS5Cl", "Na3SbS4", "Na3BrO", "NaAl11O17", "Li7P3S11", "Li2OHCl", "Na4Zr2Si3O12"])
# "NaSICON", , "LLZTO"

subscript_electrolyte_labels = []
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")  

for electrolyte in electrolyte_labels:
    subscripted_compound = electrolyte.translate(SUB)
    subscript_electrolyte_labels.append(subscripted_compound)

fontsize_ = 12
plt.rc('font', size=fontsize_) 


vol_min = np.min(vol_vals)
vol_max = np.max(vol_vals)

stress_min = np.min(stress_vals)
stress_max = np.max(stress_vals)

energy_min = np.min(energy_vals)
energy_max = np.max(energy_vals)


if (reservoir_used):
    if (open_metal_vol_included):
        cumul_ax[0].set_ylim(top=0.5)

        if ((vol_min < 0) and (vol_max > 0)): cumul_ax[1].set_ylim(bottom=vol_min*1.4, top=vol_max*2)
        elif (vol_min > 0): cumul_ax[1].set_ylim(bottom=0, top=vol_max*1.2)
        elif (vol_max < 0): cumul_ax[1].set_ylim(bottom=vol_min*1.2, top=0)

        if ((stress_min < 0) and (stress_max > 0)): 
            print(f"stress vals pos and neg")
            cumul_ax[2].set_ylim(bottom=stress_min*1.4, top=stress_max*1.4)

        elif (stress_min > 0): 
            print(f"stress vals pos: {stress_max*1.2}")
            cumul_ax[2].set_ylim(0,top=stress_max*1.2)

        elif (stress_max < 0): 
            print(f"stress vals neg")
            cumul_ax[2].set_ylim(bottom=stress_min*1.2)

    else:
        cumul_ax[0].set_ylim(top=0.5)

        if ((vol_min < 0) and (vol_max > 0)): cumul_ax[1].set_ylim(bottom=vol_min*1.4, top=vol_max*2)
        elif (vol_min > 0): cumul_ax[1].set_ylim(bottom=0, top=vol_max*1.2)
        elif (vol_max < 0): cumul_ax[1].set_ylim(bottom=vol_min*1.2, top=0)

        if ((stress_min < 0) and (stress_max > 0)): 
            print(f"stress vals pos and neg")
            cumul_ax[2].set_ylim(bottom=stress_min*1.4, top=stress_max*1.4)

        elif (stress_min > 0): 
            print(f"stress vals pos: {stress_max*1.2}")
            cumul_ax[2].set_ylim(0,top=stress_max*1.2)

        elif (stress_max < 0): 
            print(f"stress vals neg")
            cumul_ax[2].set_ylim(bottom=stress_min*1.2)

else:
    cumul_ax[0].set_ylim(top=0.01, bottom= -0.095)

figure = plt.gcf() # get current figure
figure.set_size_inches(14, 16)

cumul_ax[0].set_ylabel("Reaction \nEnergy \n(eV/atom)", size = fontsize_)

if (open_metal_vol_included): cumul_ax[1].set_ylabel("Volume change \nrelative to \nreactants (%)", size = fontsize_)
else: cumul_ax[1].set_ylabel("Volume change \nrelative to \nSE (%)", size = fontsize_)

cumul_ax[2].set_ylabel(r"$ε_{trans}E_{youngs}$ (GPa)", size = fontsize_)

cumul_ax[0].set_xticks([])
cumul_ax[0].set_xticklabels([])
cumul_ax[1].set_xticks([])
cumul_ax[1].set_xticklabels([])
# cumul_ax[2].set_xticks(np.arange(0, len(subscript_electrolyte_labels), 1))
# cumul_ax[2].set_xticklabels(subscript_electrolyte_labels, size = fontsize_, rotation=80)
print(f"cumul_ax_x_ticks_labels: {cumul_ax_x_ticks_labels}")
print(f"cumul_ax_x_ticks: {cumul_ax_x_ticks}")
cumul_ax[2].set_xticks(cumul_ax_x_ticks)
cumul_ax[2].set_xticklabels(cumul_ax_x_ticks_labels, size = fontsize_, rotation=80)

if (reservoir_used):
    if (open_metal_vol_included):
        cumul_ax[0].set_yticks(np.arange(-4, 1, 1))
        cumul_ax[0].set_yticklabels(["-4", "-3", "-2", "-1", "0"], size = fontsize_)
        cumul_ax[1].set_yticks(np.arange(-50, 11, 10))
        cumul_ax[1].set_yticklabels(["", "-40", "", "-20", "", "0", ""], size = fontsize_)
        cumul_ax[2].set_yticks(np.arange(-125, 101, 25))
        cumul_ax[2].set_yticklabels(["", "-100", "", "-50", "", "0", "", "50", "", "100"], size = fontsize_)

    else:
        cumul_ax[0].set_yticks(np.arange(-4, 1, 1))
        cumul_ax[0].set_yticklabels(["-4", "-3", "-2", "-1", "0"], size = fontsize_)
        cumul_ax[1].set_yticks(np.arange(0, 151, 25))
        cumul_ax[1].set_yticklabels(["0", "", "50", "", "100", "", "150"], size = fontsize_)
        cumul_ax[2].set_yticks(np.arange(0, 151, 25))
        cumul_ax[2].set_yticklabels(["0", "", "50", "", "100",  "", "150"], size = fontsize_)

else: 
    cumul_ax[0].set_yticks(np.arange(-0.1, 0.02, 0.02))
    cumul_ax[0].set_yticklabels(["-0.1", "-0.08", "-0.06", "-0.04", "-0.02", "0"], size = fontsize_)
    #cumul_ax[1].set_yticks(np.arange(0, 126, 25))
    #cumul_ax[1].set_yticklabels(["0", "", "50", "", "100", ""], size = fontsize_)
    cumul_ax[2].set_yticks(np.arange(0, 251, 50))
    cumul_ax[2].set_yticklabels(["0", "50", "100", "150", "200", "250"], size = fontsize_)

cumul_ax[0].set_xlim(-1,len(subscript_electrolyte_labels))
cumul_ax[1].set_xlim(-1,len(subscript_electrolyte_labels))
cumul_ax[2].set_xlim(-1,len(subscript_electrolyte_labels))

for subax in ax:
    subax.tick_params(labelsize = fontsize_)

circ1 = mpatches.Patch( facecolor="white",alpha=0.7,hatch="////",label='-1V')
circ2 = mpatches.Patch( facecolor="white",alpha=0.7,hatch=None,label='0V')
circ3 = mpatches.Patch(facecolor="white",alpha=0.7,hatch="\\\\\\\\",label='1V')
circ4 = mpatches.Patch( facecolor=blue_scalarMap.to_rgba(2),alpha=1,label='Sulfides')
circ5 = mpatches.Patch( facecolor=orange_scalarMap.to_rgba(2),alpha=1,label='Oxides')
circ6 = mpatches.Patch(facecolor=green_scalarMap.to_rgba(2),alpha=1,label='Halides')

#hatch_options = ["////", None, "\\\\\\\\"]

#cumul_ax[0].legend(handles = [circ1, circ2, circ3],loc=4, fontsize = "xx-small")
#cumul_ax[1].legend(handles = [circ1, circ2, circ3],loc=2, fontsize = "xx-small")

cumul_fig.savefig(f"{plot_dirpath}/volumechange_noLi_init.png",dpi=300)
print(plot_handles_cumul)

cumul_ax[2].legend(handles = [circ1, circ2, circ3, circ4, circ5, circ6],loc=2, fontsize = "xx-small")

for i in range(len(plot_handles_cumul_allvoltages)):
    legend = cumul_ax[1].legend(plot_handles_cumul_allvoltages[i].values(), plot_handles_cumul_allvoltages[i].keys(), bbox_to_anchor=(0,0), loc="center", fontsize = "xx-small")
    legend.get_frame().set_alpha(1)
    cumul_fig.savefig(f"{plot_dirpath}/volumechange_noLi_init_legend_{desired_voltages[i]}V.png",dpi=300)


plt.rc('font', size=fontsize_) 

print(f"trans_strains: {trans_strains}")
for i in range(len(trans_strains)):
    cumul_line_ax.plot(reaction_voltages, trans_strains[i], color = cumulativeplots_scalarMap.to_rgba(i))

print(trans_strains)
cumul_line_ax.set_ylabel(r"$ε_{trans}$", size = fontsize_)
cumul_line_ax.set_title("Transformation strains of all materials", size = fontsize_)
cumul_line_ax.set_xlim(np.min(reaction_voltages)-0.25,np.max(reaction_voltages)+0.25)
cumul_line_ax.set_xlabel("Voltage of Li or Na relative to Li/Li+ or Na/Na+ (V)", size = fontsize_)
cumul_line_fig.savefig(f"{plot_dirpath}/trans_strain_allvoltages.png",dpi=300)
legend = cumul_line_fig.legend(plot_line_handles.values(), plot_line_handles.keys(), fontsize = "xx-small")
legend.get_frame().set_alpha(1)
cumul_line_fig.savefig(f"{plot_dirpath}/trans_strain_allvoltages_legend.png",dpi=300)
cumul_line_ax.clear()

print(f"vol_ratios: {vol_ratios}")
for i in range(len(vol_ratios)):
    cumul_line_ax.plot(reaction_voltages, vol_ratios[i], color = cumulativeplots_scalarMap.to_rgba(i))
    cumul_line_ax.scatter(reaction_voltages, vol_ratios[i], color = cumulativeplots_scalarMap.to_rgba(i), alpha = 1, label=f"{labelval}", s=8, marker='s')

print(vol_ratios)
cumul_line_ax.set_ylabel(r"$V_{f} / V_{i}$", size = fontsize_)
cumul_line_ax.set_title("Final vs. initial volume ratio of all materials", size = fontsize_)
cumul_line_ax.set_xlim(np.min(reaction_voltages)-0.25,np.max(reaction_voltages)+0.25)
cumul_line_ax.set_xlabel("Voltage of Li or Na relative to Li/Li+ or Na/Na+ (V)", size = fontsize_)
cumul_line_fig.savefig(f"{plot_dirpath}/vol_ratios_allvoltages.png",dpi=300)
legend = cumul_line_fig.legend(plot_line_handles.values(), plot_line_handles.keys(), fontsize = "xx-small")
legend.get_frame().set_alpha(1)
cumul_line_fig.savefig(f"{plot_dirpath}/vol_ratios_allvoltages_legend.png",dpi=300)

fontsize_ = 20
print(f"reaction_plot_handles: {reaction_plot_handles}")
reactions_ax[0].set_ylabel(r"Reaction energy (eV/atom)", size = fontsize_)
reactions_ax[1].set_ylabel(r"$V_{f} / V_{i}$ (%)", size = fontsize_)
reactions_ax[0].set_title(r"$Li_{7}P_{3}S_{11}$ reaction energies and \nvolume changes with lithium metal at various voltages", size = fontsize_)
reactions_ax[0].set_xlim(np.min(reaction_voltages)-0.25,np.max(reaction_voltages)+0.25)
reactions_ax[1].set_xlim(np.min(reaction_voltages)-0.25,np.max(reaction_voltages)+0.25)
reactions_ax[1].set_xlabel("Voltage of Li relative to Li/Li+ (V)", size = fontsize_)
reactions_fig.savefig(f"{plot_dirpath}/reaction_shifted.png",dpi=300)
legend = reactions_fig.legend(reaction_plot_handles.values(), reaction_plot_handles.keys(), fontsize = "small")
legend.get_frame().set_alpha(1)
reactions_fig.savefig(f"{plot_dirpath}/reaction_shifted_legend.png",dpi=300)


workbook.close()