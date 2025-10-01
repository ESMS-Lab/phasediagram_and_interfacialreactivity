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
#from mendeleev import *


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
            coeff = product_split[0]
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
        
        if (reactant == reactant_side_comp_coeff[-1]): new_reactant_str += " -> "
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


def assign_vol_changes(volfilepath, reactfilepath):
        
    volfile = open("volume_changes_stable_products/" + volfilepath, "r")
    #volfile = open("test_volume_changes/" + volfilepath, "r")
    volfile_lines = [line for line in volfile]
    
    #reactfile = open("all_reactions/" + reactfilepath, "r")
    #reactfile = open("test_reactions/" + reactfilepath, "r")
    reactfile = open("no_baremetal_reactions/" + reactfilepath, "r")
    reactfile_lines = [line for line in reactfile]

    print(f"reactfile: {reactfilepath}")
    
    reactfile_split = (reactfilepath.split("_"))
    elems = reactfile_split[:-2]
    elems_joined =  "_".join(elems)
    voltage = reactfile_split[-2]
    #vol_pot_name = "volumes_by_voltage/" + elems_joined
    vol_pot_name = "test_volumes_by_voltage/" + elems_joined
    vol_pot_name += "_volumes_by_voltage.txt"
    vol_by_pot_file = open(vol_pot_name, 'a')
    
    compound_dict = parse_volfile(volfile_lines)
        
    vol_by_pot_file.write(f"elements: {elems_joined} at {voltage}V\n") 
    vol_by_pot_file.write(f"--------------------------------------\n") 

    reac_equ_dict = {}
    reac_energy_dict = {}
    reac_electrolyte_mol_dict = {}
    
    for react_line in reactfile_lines:
        print(f"react_line: {react_line}")
        if ("--------" in react_line) or ("decomposition reactions" in react_line) or react_line.isspace(): pass
        else:
            #print(f"react_line: {react_line}")
            react_dict, prod_dict, react_energy = parse_reaction(react_line)

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

            ### generating most stable version of compound ###
            '''
            reactants = []
            products = []
            print(f"len(react_dict): {len(react_dict)} len(prod_dict): {len(prod_dict)}")
            
            for compound, coeff in react_dict.items():
                reactants.append(get_most_stable_entry(compound, all_elems))
            
            for compound, coeff in prod_dict.items():
                products.append(get_most_stable_entry(compound, all_elems))
            '''
            #print(f"reactants: {reactants} products: {products}")
            

            reactant_vol = 0
            electrolyte_coeff = 0
            
            for compound, coeff in react_dict.items():
                print(f"reactant compound: {compound}")
                if (compound not in compound_dict): KeyError("Key 'compound' is not in dictionary 'compound_dict'")
                elif ((compound == "Li1") or (compound == "Na1")): 
                     pass # include this if statement in case of omitting volume of open element
                else: 
                    molar_volume = compound_dict[compound]
                    reactant_vol += float(molar_volume) * float(coeff)
                                
                compound_split = seperate_string_number(compound)
                compound_elements = [str(elem) for elem in compound_split if not elem.isdigit()]
                #print(f"compound_elements: {compound_elements}")
                #print(f"elems: {elems}")
                if (set(compound_elements) == set(elems)): electrolyte_coeff = coeff

            #print(f"end of react")
            product_vol = 0
            for compound, coeff in prod_dict.items():
                #print(f"product compound: {compound}")
                if (compound not in compound_dict): KeyError("Key 'compound' is not in dictionary 'compound_dict'")
                else: 
                    molar_volume = compound_dict[compound]
                    #print(f"product compound: {compound}")
                    #print(f"product molar_volume: {molar_volume}   product coeff: {coeff}")
                    product_vol += float(molar_volume) * float(coeff)

            volume_change = product_vol - reactant_vol
            print(f"reactant_vol: {reactant_vol}   product_vol: {product_vol}   volume_change: {volume_change}")
            rel_vol_change_percentage = ( volume_change / reactant_vol ) * 100
            
            react_line_strip = " ".join(react_line.split()[:-2])
            reac_equ_dict[react_line_strip] = rel_vol_change_percentage
            reac_energy_dict[react_line_strip] = react_energy
            reac_electrolyte_mol_dict[react_line_strip] = electrolyte_coeff
            vol_by_pot_file.write(react_line_strip + f" volume change: {round(volume_change, 4)}, relative vol change: {round(rel_vol_change_percentage, 4)}% \n ") 

    vol_by_pot_file.write(f"\n\n") 

    return reac_equ_dict, reac_energy_dict, reac_electrolyte_mol_dict


def write_to_spreadsheet(reactions, vol_changes, voltage, electrolyte_name):
    # Workbook() takes one, non-optional, argument 
    # which is the filename that we want to create.
    workbook = xlsxwriter.Workbook(f"{electrolyte_name}_decomp_reactions.xlsx")
    row = 1 + voltage_shift 
    
    for k in range(len(reacts_sorted)):
        col = 1
        workbook.write(row, col, vol_changes_sorted[k])
     

def find_max_min_in_dict(dictionary, youngs_mod=[]):
    val_max = 0
    val_min = 0
    ### setting y bounds for all volume plots ###
    for i in range(len(dictionary)):
        for j in range(len(dictionary[0])):
            print(f"i: {i} j: {j}")
            for reaction, val_change in dictionary[i][j].items():
                if (len(youngs_mod) == 0): pass
                else: val_change = np.cbrt(val_change) * youngs_mod[i]

                if (val_change > val_max): val_max = val_change
                elif (val_change < val_min): val_min = val_change

    return val_min, val_max


#dirpath = os.getcwd() + "/all_reactions/"
#dirpath = os.getcwd() + "/test_reactions/"
dirpath = os.getcwd() + "/no_baremetal_reactions/"
#plot_dirpath = "test_plots"
#plot_dirpath = "plots"
plot_dirpath = "no_baremetal_plots"

reaction_elems = np.asarray([set(["Li", "Ge", "P", "S"]), set(["Li", "In", "Cl"]), set(["Li", "La", "Zr", "O"]), 
    set(["Li", "La", "Zr", "Ta", "O"]), set(["Li", "P", "S", "Cl"]), set(["Na", "Sb", "S"]), set(["Na", "Zr", "Si", "P", "O"]), 
    set(["Na", "Br", "O"]), set(["Al", "Na", "O"]), set(["Li", "P", "S"]), set(["Li", "O", "H", "Cl"])])

electrolyte_compounds = np.asarray(["Li10GeP2S12", "LiInCl6", "Li7La3Zr2O12", "Li13O24La6Zr3Ta",  
   "Li6PS5Cl", "Na3SbS4", "Na3Zr2Si2PO12", "Na3BrO", "NaAl11O17", "Li7P3S11", "Li2OHCl"])


reaction_dicts = np.empty((len(electrolyte_compounds),21),dtype=object)
reaction_energies = np.empty((len(electrolyte_compounds),21),dtype=object)
reaction_elec_mols = np.empty((len(electrolyte_compounds),21),dtype=object)

electrolyte_youngs_mod = np.asarray([37.19, 100, 149.8, 100, 22 , 33.9, 56, 100, 100, 100, 100])
#reaction_voltages = np.round(np.arange(-0.5, 0.51, 0.05), 2)
reaction_voltages = np.round(np.arange(-1.0, 1.01, 0.1), 2)

for file in os.listdir(dirpath):

    input_file_split = (file.split("_"))
    input_elems = input_file_split[:-2]
    # sub_dirpath = os.getcwd() + "/volume_changes_stable_products/"
    sub_dirpath = os.getcwd() + "/test_volume_changes/"

    for sub_file in os.listdir(sub_dirpath):

        print(f"input_elems: {input_elems}")
        print(f"reaction_elems: {reaction_elems}")
        # if (set(input_elems) in reaction_elems):
        sub_file_split = (sub_file.split("_"))
        sub_elems = sub_file_split[:-2]

        if (set(sub_elems) == set(input_elems)):
            
            print(f"sub_elems: {sub_elems}    input_elems: {input_elems}")
            reac_eqn_dict, react_energy_dict, react_elec_mol_dict = assign_vol_changes(sub_file, file)
            filename_split = file.split("_")
            elems = set(filename_split[:-2])
            
            voltage = float(filename_split[-2].strip("volts"))
            print(f"voltage: {voltage}")
            elem_i = np.where(reaction_elems == elems)[0][0]
            voltage_i = np.where(reaction_voltages == voltage)[0][0]
            reaction_dicts[elem_i][voltage_i] = reac_eqn_dict
            reaction_energies[elem_i][voltage_i] = react_energy_dict
            reaction_elec_mols[elem_i][voltage_i] = react_elec_mol_dict


print(f"reaction_dicts: {reaction_dicts}")

### setting y bounds for all volume plots ###
vol_min, vol_max = find_max_min_in_dict(reaction_dicts)

### setting y bounds for stress plots ###
stress_min, stress_max = find_max_min_in_dict(reaction_dicts, youngs_mod=electrolyte_youngs_mod)
energy_min, energy_max = find_max_min_in_dict(reaction_energies)


electrolyte_labels = np.asarray(["LGPS", "Li3InCl6", "LLZO", "LLZTO",  
   "LPSCl", "Na3SbS4", "NaSICON", "Na3BrO", "NaAl11O17", "Li7P3S11", "Li2OHCl"])

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
cmap_name = "gnuplot2"
x_pos = 0
cumulative_plot_handles = {}
start_val = 0
stop_val = len(reaction_dicts) + 1
cmap_in = plt.get_cmap(cmap_name)
norm_in = mpl.colors.Normalize(vmin=start_val-1, vmax=stop_val+1)
cumulativeplots_scalarMap = cm.ScalarMappable(norm=norm_in, cmap=cmap_in)
cumul_fig, cumul_ax = plt.subplots((3), figsize=(6,6), constrained_layout = True)
cumul_line_fig, cumul_line_ax = plt.subplots(figsize=(8,6), constrained_layout = True)
trans_strains = np.zeros((len(reaction_dicts), len(reaction_dicts[0])))

# Creating figure objects for plots for each material
fig, ax = plt.subplots((3), figsize=(4,6), constrained_layout = True)
fontsize_ = 34
plot_handles_cumul = {}
plot_line_handles = {}

### voltage for plot of transformation stran, volume change, and reaction energy ###
### spanning all compounds ###
desired_voltage = 0
j_cumulative = np.where(reaction_voltages == desired_voltage)[0][0]


### plotting data ###
for i in range(len(reaction_dicts)):
    
    print(reaction_dicts.size)
    print(f"\n\ni: {i}  electrolyte_compounds[i]: {electrolyte_compounds[i]}")
    
    ### creating list (non-repetitive) of all reactions for given solid ###
    ### electrolyte, removing 0-volume change reactions ###
    ax_reactions = []
    ax_energies = np.zeros((len(reaction_dicts[0])))
    ax_electrolyte = " ".join(reaction_elems[i])
    plot_handles = {} #Uncomment for multiple plots

    for j in range(len(reaction_dicts[0])):
        # print(f"j: {j}")
        for reaction, volume_change in reaction_dicts[i][j].items():
                # print(f"reaction: {reaction}  volume_change: {volume_change}   reaction_energies[i][j][str(reaction)]: {reaction_energies[i][j][str(reaction)]}")
                                
                if (volume_change == 0) or (reaction_elec_mols[i][j][str(reaction)] < 1):
                    # print(f"pass")
                    pass

                else: 
                    # print(f"ACCEPTED reaction: {reaction}  volume_change: {volume_change}   reaction_energies[i][j][str(reaction)]: {reaction_energies[i][j][str(reaction)]}")
                    ax_reactions.append(str(reaction))
                    ax_energies[j] = (reaction_energies[i][j][str(reaction)])
    
    ### removing duplicate reactions with permutations of compounds  ###
    set_reacts = []
    noduplicates_ax_reactions = []
    ax_reactions = np.asarray(list(set(ax_reactions)))
    subscripted_ax_reactions = []

    for j in range(len(ax_reactions)):
        reaction = ax_reactions[j] 
        #energy = ax_energies[j] 
        # print(f"reaction idx: {reaction}")
        reactant_dict, product_dict, react_energy = parse_reaction(reaction)
        # is_redox = check_if_redox(reactant_dict, product_dict)
        new_reac_set = set()
        
        # print(f"set reaction: {reaction} ")

        for compound, coeff in product_dict.items():
            new_reac_set.add(compound)
        
        if new_reac_set not in set_reacts: 
            # print(f"ACCEPTED set reaction: {reaction} ")
            set_reacts.append(new_reac_set)
            noduplicates_ax_reactions.append(reaction)

    print(f"set_reacts:")
    print(set_reacts)
    
    ax_reactions = np.asarray(noduplicates_ax_reactions)

    for react in ax_reactions:        
        subscripted_ax_reactions.append(make_reaction_subscripted(react)) 
    
    # uncomment for multiple plots
    start_val = 0
    stop_val = len(subscripted_ax_reactions)
    cmap_in = plt.get_cmap(cmap_name)
    norm_in = mpl.colors.Normalize(vmin=start_val-1, vmax=stop_val+1)
    scalarMap = cm.ScalarMappable(norm=norm_in, cmap=cmap_in)

    ### creating bars ###
    for j in range(len(reaction_dicts[0])):  ### TAB BACK THE REST OF THIS FOR LOOP
        
        ax_voltage = reaction_voltages[j] 

        '''
        start_val = 0
        stop_val = len(ax_reactions) + 1
        cmap_in = plt.get_cmap(cmap_name)
        norm_in = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
        scalarMap = cm.ScalarMappable(norm=norm_in, cmap=cmap_in)
        '''
        
        color_i = 0
        
        vol_changes_unsorted = []
        reacts_unsorted = []
        energies_unsorted = []
        elec_mols_unsorted = []
        
        ### removing 0-volume reactions and selecting min-energy reaction ###
        ### at each given voltage (multiple can occur) ###
        if (len(reaction_dicts[i][j]) > 1):
            min_energy = sys.maxsize
            min_e_reaction = ""
            min_e_vol_change = 0
            min_e_elec_mols = 0
            
            for reaction, volume_change in reaction_dicts[i][j].items():
                # print(f"reaction: {reaction}  reaction_energies[i][j][reaction]: {reaction_energies[i][j][reaction]}")
                if (volume_change == 0) or (reaction_elec_mols[i][j][reaction] < 1): pass
                elif (reaction_energies[i][j][reaction] < min_energy):
                    min_energy = reaction_energies[i][j][reaction]
                    min_e_reaction = reaction
                    min_e_vol_change = volume_change
                    min_e_elec_mols = reaction_elec_mols[i][j][reaction]

            if not min_e_reaction: pass 
            else:
                vol_changes_unsorted.append(min_e_vol_change)
                reacts_unsorted.append(min_e_reaction)      
                energies_unsorted.append(min_energy)
                elec_mols_unsorted.append(min_e_elec_mols)

        else:
            for reaction, volume_change in reaction_dicts[i][j].items():
                if (volume_change == 0) or (reaction_elec_mols[i][j][reaction] != 1):
                    print(f"pass")
                    pass
                else: 
                    vol_changes_unsorted.append(volume_change)
                    reacts_unsorted.append(reaction)
                    energies_unsorted.append(reaction_energies[i][j][reaction])
                    elec_mols_unsorted.append(reaction_elec_mols[i][j][reaction])
                    print(f"reaction: {reaction}")

        reacts_sorted = [react for _, react in sorted(zip(vol_changes_unsorted, reacts_unsorted))]
        print(f"reacts_sorted: {reacts_sorted}")
        energies_sorted = [react for _, react in sorted(zip(vol_changes_unsorted, energies_unsorted))]
        elec_mols_sorted = [react for _, react in sorted(zip(vol_changes_unsorted, elec_mols_unsorted))]
        scaled_volumes_sorted = np.sort(vol_changes_unsorted) 


        vol_changes_unsorted = []
        reacts_unsorted = []
        energies_unsorted = []
        elec_mols_unsorted = []
        
        ### removing 0-volume reactions and selecting min-energy reaction ###
        ### at each given voltage (multiple can occur) ###
        if (len(reaction_dicts[i][j_cumulative]) > 1):
            min_energy = sys.maxsize
            min_e_reaction = ""
            min_e_vol_change = 0
            min_e_elec_mols = 0
            
            for reaction, volume_change in reaction_dicts[i][j_cumulative].items():
                # print(f"reaction: {reaction}  reaction_energies[i][j][reaction]: {reaction_energies[i][j][reaction]}")
                if (volume_change == 0) or (reaction_elec_mols[i][j_cumulative][reaction] < 1): pass
                elif (reaction_energies[i][j_cumulative][reaction] < min_energy):
                    min_energy = reaction_energies[i][j_cumulative][reaction]
                    min_e_reaction = reaction
                    min_e_vol_change = volume_change
                    min_e_elec_mols = reaction_elec_mols[i][j_cumulative][reaction]

            if not min_e_reaction: pass 
            else:
                vol_changes_unsorted.append(min_e_vol_change)
                reacts_unsorted.append(min_e_reaction)      
                energies_unsorted.append(min_energy)
                elec_mols_unsorted.append(min_e_elec_mols)

        else:
            for reaction, volume_change in reaction_dicts[i][j_cumulative].items():
                if (volume_change == 0) or (reaction_elec_mols[i][j_cumulative][reaction] != 1):
                    print(f"pass")
                    pass
                else: 
                    vol_changes_unsorted.append(volume_change)
                    reacts_unsorted.append(reaction)
                    energies_unsorted.append(reaction_energies[i][j][reaction])
                    elec_mols_unsorted.append(reaction_elec_mols[i][j][reaction])
                    print(f"reaction: {reaction}")

        reacts_sorted_cumul = [react for _, react in sorted(zip(vol_changes_unsorted, reacts_unsorted))]
        print(f"reacts_sorted: {reacts_sorted}")
        energies_sorted_cumul = [react for _, react in sorted(zip(vol_changes_unsorted, energies_unsorted))]
        elec_mols_sorted_cumul = [react for _, react in sorted(zip(vol_changes_unsorted, elec_mols_unsorted))]
        scaled_volumes_sorted_cumul = np.sort(vol_changes_unsorted) 

        ### (deprecated) if more than one reaction at given voltage, scaling ###
        ###   volume changes according to abundance based on energies of reactions ###
        T = 300
        kB = 8.6173e-5
        total_vol_change = 0
        total_stress = 0


        ### plotting bars according to volume change for each reaction ###
        if (not reacts_sorted):
            labelval = f"None"
            bar_handle = ax[0].bar(x_pos, 0, width = 1, color = scalarMap.to_rgba(j), alpha = 1, label=f"{labelval}")
            bar_handle = cumul_ax[0].bar(x_pos, 0, width = 1, color = scalarMap.to_rgba(i), alpha = 1, label=f"{labelval}")
            print(f"react sorted empty")
            
        for k in range(len(reacts_sorted)):
            
            reactant_dict, product_dict, react_energy = parse_reaction(reacts_sorted[k])
            new_reac_set = set()
            for compound, coeff in product_dict.items():
                new_reac_set.add(compound)

            color_i = np.where(np.asarray(set_reacts) == new_reac_set)[0][0]
            
            if (scaled_volumes_sorted[k] == 0): 
                print(f"error: vol_changes_sorted[k] should have no 0 elements")
                pass
            else:
                
                labelval = f"{subscripted_ax_reactions[color_i]}"
                bar_handle = ax[1].bar(ax_voltage, scaled_volumes_sorted[k], width = round(reaction_voltages[1] - reaction_voltages[0], 2), color = scalarMap.to_rgba(color_i), alpha = 1, label=f"{labelval}")
                 
                if (bar_handle not in plot_handles): plot_handles[labelval] = bar_handle
                
                total_vol_change += scaled_volumes_sorted[k]
                bar_handle = ax[0].bar(ax_voltage, np.round(ax_energies[j],3), width = round(reaction_voltages[1] - reaction_voltages[0], 2), color = scalarMap.to_rgba(color_i), alpha = 1, label=f"{labelval}")
                bar_line_handle_cumul = cumul_line_ax.scatter(ax_voltage, np.cbrt(scaled_volumes_sorted[k] / 100), color = cumulativeplots_scalarMap.to_rgba(i), alpha = 1, label=f"{labelval}", s=8, marker='s')
                trans_strains[i][j] = np.cbrt(scaled_volumes_sorted[k] / 100)
                if (bar_line_handle_cumul not in plot_line_handles): plot_line_handles[subscript_electrolyte_compounds[i]] = bar_line_handle_cumul

        for k in range(len(reacts_sorted_cumul)):
            
            reactant_dict, product_dict, react_energy = parse_reaction(reacts_sorted_cumul[k])
            new_reac_set = set()
            for compound, coeff in product_dict.items():
                new_reac_set.add(compound)

            color_i = np.where(np.asarray(set_reacts) == new_reac_set)[0][0]
            
            if (scaled_volumes_sorted_cumul[k] == 0): 
                print(f"error: vol_changes_sorted[k] should have no 0 elements")
                pass
            else:
                print(f"bar for i: {i}")
                print(f"x_pos: {x_pos}")
                print(f"scaled_volumes_sorted_cumul[k]: {scaled_volumes_sorted_cumul[k]}")
                labelval = f"{subscripted_ax_reactions[color_i]}"
                bar_handle_cumul = cumul_ax[1].bar(x_pos, scaled_volumes_sorted_cumul[k], width = 1, color = cumulativeplots_scalarMap.to_rgba(i), alpha = 1, label=f"{labelval}")

                if (bar_handle_cumul not in plot_handles_cumul): plot_handles_cumul[labelval] = bar_handle_cumul

                total_vol_change += scaled_volumes_sorted_cumul[k]
                bar_handle_cumul = cumul_ax[0].bar(x_pos, np.round(ax_energies[j_cumulative],3), width = 1, color = cumulativeplots_scalarMap.to_rgba(i), alpha = 1, label=f"{labelval}")
                
        for k in range(len(reacts_sorted)):

            reactant_dict, product_dict, react_energy = parse_reaction(reacts_sorted[k])
            new_reac_set = set()
            for compound, coeff in product_dict.items():
                new_reac_set.add(compound)

            color_i = np.where(np.asarray(set_reacts) == new_reac_set)[0][0]

            if (scaled_volumes_sorted[k] == 0): 
                print(f"error: vol_changes_sorted[k] should have no 0 elements")
                pass
            else:
                ###  calculating change in effective fracture toughness  ###
                ### redox active ==> reaction at tip, chemical ==> reaction along crack length ###
                labelval = f"{subscripted_ax_reactions[color_i]}"
                #alpha = 0.42
                #thickness = 3.75e-8
                deltaK_div_thickness =  np.cbrt(scaled_volumes_sorted[k]) * electrolyte_youngs_mod[i]
                bar_handle = ax[2].bar(ax_voltage, deltaK_div_thickness, width = round(reaction_voltages[1] - reaction_voltages[0], 2), color = scalarMap.to_rgba(color_i), alpha = 1, label=f"{ax_reactions[color_i]}")
                
        for k in range(len(reacts_sorted_cumul)):

            reactant_dict, product_dict, react_energy = parse_reaction(reacts_sorted_cumul[k])
            new_reac_set = set()
            for compound, coeff in product_dict.items():
                new_reac_set.add(compound)

            color_i = np.where(np.asarray(set_reacts) == new_reac_set)[0][0]

            if (scaled_volumes_sorted_cumul[k] == 0): 
                print(f"error: vol_changes_sorted[k] should have no 0 elements")
                pass
            else:
                ###  calculating change in effective fracture toughness  ###
                ### redox active ==> reaction at tip, chemical ==> reaction along crack length ###
                labelval = f"{subscripted_ax_reactions[color_i]}"
                #alpha = 0.42
                #thickness = 3.75e-8
                deltaK_div_thickness =  np.cbrt(scaled_volumes_sorted_cumul[k]) * electrolyte_youngs_mod[i]
                bar_handle_cumul = cumul_ax[2].bar(x_pos, deltaK_div_thickness, width = 1, color = cumulativeplots_scalarMap.to_rgba(i), alpha = 1, label=f"{ax_reactions[color_i]}")
                

        ### setting yscale for axis object ###
        if (total_vol_change > vol_max): vol_max = total_vol_change
        elif (total_vol_change < vol_min): vol_min = total_vol_change

        if (total_stress > stress_max): stress_max = total_stress   
        elif (total_stress < stress_min): stress_min = total_stress
    
    x_pos += 1
        
    print("\n\n")
    fontsize_ = 12
    plt.rc('font', size=fontsize_) 

    ### Plotting reactions for indivdual materials ###
    print(f"i: {i}   ax_energies: {ax_energies}, len(ax_energies): {len(ax_energies)}")
    print(f'vol_min: {vol_min}')
    print(f'vol_max: {vol_max}')
    
    ax[0].set_title(subscript_electrolyte_labels[i])
    ax[0].set_xlim(np.min(reaction_voltages)-0.25,np.max(reaction_voltages)+0.25)
    ax[0].set_ylabel("Reaction \nEnergy \n(eV/atom)", size = fontsize_)

    ax[1].set_xlim(np.min(reaction_voltages)-0.25,np.max(reaction_voltages)+0.25)
    ax[1].set_ylabel("Relative \nvolume \nchange (%)", size = fontsize_)
    
    ax[2].set_xlim(np.min(reaction_voltages)-0.25,np.max(reaction_voltages)+0.25)
    ax[2].set_ylabel(r"$ε_{trans}E_{youngs}$", size = fontsize_)
    ax[2].set_xlabel("Potential of Li or Na (eV)")
    
    figure = plt.gcf() # get current figure
    figure.set_size_inches(10, 6)

    fig.savefig(f"{plot_dirpath}/stress_electrolyte_{electrolyte_compounds[i]}.png",dpi=300)

    legend = ax[1].legend(plot_handles.values(), plot_handles.keys(), bbox_to_anchor=(-0.15, 0), loc="lower left", fontsize = "x-small")
    legend.get_frame().set_alpha(1)            
    figure = plt.gcf() # get current figure
    fig.savefig(f"{plot_dirpath}/legend_volchange_electrolyte_{electrolyte_compounds[i]}.png",dpi=300)
    
    ax[0].clear()
    ax[1].clear()
    ax[2].clear()

### Plotting cumulative plots ###
electrolyte_labels = np.asarray(["LGPS", "Li3InCl6", "LLZO", "LLZTO",  
   "LPSCl", "Na3SbS4", "NaSICON", "Na3BrO", "NaAl11O17", "Li7P3S11", "Li2OHCl"])

subscript_electrolyte_labels = []
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")  

for electrolyte in electrolyte_labels:
        subscripted_compound = electrolyte.translate(SUB)
        subscript_electrolyte_labels.append(subscripted_compound)

fontsize_ = 12
plt.rc('font', size=fontsize_) 

cumul_ax[0].set_ylim(top=0.5)
figure = plt.gcf() # get current figure
figure.set_size_inches(14, 16)
cumul_ax[1].set_ylabel("Relative \nvolume \nchange (%)", size = fontsize_)
cumul_ax[2].set_ylabel(r"$ε_{trans}E_{youngs}$", size = fontsize_)
cumul_ax[0].set_ylabel("Reaction \nEnergy \n(eV/atom)", size = fontsize_)

cumul_ax[0].set_xticks([])
cumul_ax[0].set_xticklabels([])
cumul_ax[1].set_xticks([])
cumul_ax[1].set_xticklabels([])
cumul_ax[2].set_xticks(np.arange(0, len(subscript_electrolyte_labels), 1))
cumul_ax[2].set_xticklabels(subscript_electrolyte_labels, size = fontsize_, rotation=80)
cumul_ax[0].set_xlim(-1,len(subscript_electrolyte_labels))
cumul_ax[1].set_xlim(-1,len(subscript_electrolyte_labels))
cumul_ax[2].set_xlim(-1,len(subscript_electrolyte_labels))

for subax in ax:
    subax.tick_params(labelsize = fontsize_)

cumul_fig.savefig(f"{plot_dirpath}/volumechange_{desired_voltage}v_noLi_init.png",dpi=300)
print(plot_handles_cumul)
legend = cumul_ax[1].legend(plot_handles_cumul.values(), plot_handles_cumul.keys(), bbox_to_anchor=(0,0), loc="center", fontsize = "xx-small")
legend.get_frame().set_alpha(1)
cumul_fig.savefig(f"{plot_dirpath}/volumechange_{desired_voltage}v_noLi_init_legend.png",dpi=300)


fontsize_ = 12
plt.rc('font', size=fontsize_) 
#cumul_line_ax[2].set_ylim(top=0.5)
#figure = plt.gcf() # get current figure
#figure.set_size_inches(14, 16)
for i in range(len(trans_strains)):
    cumul_line_ax.plot(reaction_voltages, trans_strains[i], color = cumulativeplots_scalarMap.to_rgba(i))

print(trans_strains)
cumul_line_ax.set_ylabel(r"$ε_{trans}$", size = fontsize_)
cumul_line_ax.set_title("Transformation strains of all materials", size = fontsize_)
cumul_line_ax.set_xlim(np.min(reaction_voltages)-0.25,np.max(reaction_voltages)+0.25)
cumul_line_ax.set_xlabel("Potential of Li or Na (eV)", size = fontsize_)
cumul_line_fig.savefig(f"{plot_dirpath}/trans_strain_allvoltages.png",dpi=300)
legend = cumul_line_fig.legend(plot_line_handles.values(), plot_line_handles.keys(), fontsize = "xx-small")
legend.get_frame().set_alpha(1)
cumul_line_fig.savefig(f"{plot_dirpath}/trans_strain_allvoltages_legend.png",dpi=300)
