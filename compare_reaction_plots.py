from pymatgen.ext.matproj import MPRester 
from pymatgen.analysis.phase_diagram import PhaseDiagram, GrandPotentialPhaseDiagram, PDEntry, PDPlotter
from pymatgen.analysis.interface_reactions import InterfacialReactivity, GrandPotentialInterfacialReactivity
from pymatgen.analysis.reaction_calculator import ComputedReaction
from pymatgen.core import Composition, Element
from pymatgen.io.vasp import Vasprun
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
from pymatgen.analysis.phase_diagram import PDPlotter
import numpy as np
import os
import copy
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import pandas 
from collections import OrderedDict
from itertools import combinations


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


def gather_reactions_and_e(filepath):
    filepath = os.path.join(os.getcwd(), filepath)
    reactfile = open(filepath, 'r')
    reactfile_lines = [line for line in reactfile]
    dict_string = ""
    
    for react_line in reactfile_lines:
        print(f"react_line: {react_line}")
        if ("--------" in react_line) or react_line.isspace(): pass
        elif ("bare metal" in react_line): 
            dict_string = "bare metal"
            reactant = (react_line.split())[0]
            reactant_split = seperate_string_number(reactant)
            reactant_elements = [str(elem) for elem in reactant_split if not elem.isdigit()]
            print(f"reactant_elements: {reactant_elements}")
            
        elif ("theoretical" in react_line): 
            dict_string = "theoretical"
            reactant = (react_line.split())[0]
            reactant_split = seperate_string_number(reactant)
            reactant_elements = [str(elem) for elem in reactant_split if not elem.isdigit()]
            print(f"reactant_elements: {reactant_elements}")

        else:
            react_line = react_line.split("normal")
            reaction_chunk = react_line[0]
            energy_chunk = react_line[1]
            grand_potential = (energy_chunk.split())[-1]
            react_dict, _, _ = parse_reaction(reaction_chunk)
            print(f"react_dict: {react_dict}")
            print(f"grand_potential: {grand_potential}")
            print(f"reactant_elements: {reactant_elements}")

            electrolyte_coeff = 0                
            for compound, coeff in react_dict.items():
                compound_split = seperate_string_number(compound)
                compound_elements = [str(elem) for elem in compound_split if not elem.isdigit()]
                if (set(compound_elements) == set(reactant_elements)): electrolyte_coeff = coeff 

            reaction_chunk = "".join(reaction_chunk)
            
            if np.isclose(electrolyte_coeff, 1):
                if (dict_string == "bare metal"): min_e_reac_e_baremetal = tuple([reaction_chunk, grand_potential])
                elif (dict_string == "theoretical"): min_e_reac_e_theo = tuple([reaction_chunk, grand_potential])
             
    reactfile.close()

    return min_e_reac_e_theo, min_e_reac_e_baremetal

mu_vals = np.round(np.arange(-1.0, 1.01, 0.1), 2)
unwanted_elements = []
unwanted_elements = ["Zr", "Ta", "Ge", "P", "La", "Sb", "Br", "In", "S", "P", "Al"]

outputdir = "baremetals_vs_theo_comparison"
only_exp_observed_val = False

overlap = []
only_theo = []
only_exp = []
custom_theoretical_compounds = "all_exp_elems/consolidated_theoretical_master_dupremoved_marked.txt"

electrolyte_compounds = np.asarray(["Li7La3Zr2O12", "Li10GeP2S12", "Na3SbS4", "Li3InCl6", 
   "Na3BrO", "Li6PS5Cl", "Li7P3S11", "NaAl11O17", "Li2OHCl", "Na4Zr2Si3O12"])
# "Li13O24La6Zr3Ta",  "
electrolyte_elems = np.asarray([set(["Li", "O", "La", "Zr"]), set(["Li", "Ge", "P", "S"]), set(["Na", "Sb", "S"]), set(["Li", "In", "Cl"]), 
    set(["Na", "Br", "O"]),  set(["Li", "P", "S", "Cl"]), set(["Li", "P", "S"]),  set(["Na", "Al", "O"]), set(["Li", "O", "H", "Cl"]), set(["Na", "Zr", "Si", "O"])])
#set(["Li", "O", "La", "Zr", "Ta"]),
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")  

plot_dirpath = "comparison_plots"
energy_difference = np.zeros((len(electrolyte_compounds), 2, len(mu_vals)))
reaction_array = np.empty((len(electrolyte_compounds), 2, len(mu_vals)),dtype=object)
idx = 0

print(f"energy_difference.shape: {energy_difference.shape}")

dirpath = os.getcwd()
for file in os.listdir(dirpath):
    if (".py" not in file) and ("reactions.txt" in file):
        print(f"file: {file}")
        file_split = (file.split("_"))
        elems = file_split[:-2]
        elem_idx = np.where(set(elems) == electrolyte_elems)[0][0]
        voltage = float(file_split[-2].strip("volts"))
        print(f"voltage: {voltage}")
        voltage_i = np.where(mu_vals == voltage)[0][0]
        min_e_reac_e_theo, min_e_reac_e_baremetal = gather_reactions_and_e(file)

        print(f"elem_idx: {elem_idx} voltage_i: {voltage_i}")
        energy_difference[elem_idx][0][voltage_i] = min_e_reac_e_theo[1]
        energy_difference[elem_idx][1][voltage_i] = min_e_reac_e_baremetal[1]

        reaction_array[elem_idx][0][voltage_i] = min_e_reac_e_theo[0]
        reaction_array[elem_idx][1][voltage_i] = min_e_reac_e_baremetal[0]


fig, ax = plt.subplots((1), figsize=(4,6), constrained_layout = True)
cmap_name = "gnuplot2"

for i in range(len(energy_difference)):
    print(f"electrolyte_compounds[i]: {electrolyte_compounds[i]}")
    product_combinations = []

    for j in range(len(energy_difference[0][0])):
        reaction_theo = reaction_array[i][0][j]
        reaction_baremetal = reaction_array[i][1][j]
        _, baremetal_products, _ = parse_reaction(reaction_baremetal)
        _, theo_products, _ = parse_reaction(reaction_theo)

        baremetal_reac_set = set()
        for compound in baremetal_products.items():
            baremetal_reac_set.add(compound)
        
        theo_reac_set = set()
        for compound in theo_products.items():
            theo_reac_set.add(compound)

        reac_combination = [theo_reac_set, baremetal_reac_set]
        if reac_combination not in product_combinations: product_combinations.append(reac_combination)

    product_combinations = np.asarray(product_combinations)
    start_val = -1
    stop_val = len(product_combinations)
    cmap_in = plt.get_cmap(cmap_name)
    norm_in = mpl.colors.Normalize(vmin=start_val-1, vmax=stop_val+1)
    scalarMap = cm.ScalarMappable(norm=norm_in, cmap=cmap_in)
    
    plot_handles = {}
    reactions = []
    label_vals = []

    for j in range(len(energy_difference[0][0])):
        if not np.isclose((energy_difference[i][1][j] - energy_difference[i][0][j]), 0):
            print(f"mu_vals[j]: {mu_vals[j]}")
            ### Plotting reactions for indivdual materials ###
            j_reversed = int(len(energy_difference[0][0])) - j - 1
            reaction_theo = reaction_array[i][0][j]
            reaction_baremetal = reaction_array[i][1][j]
            _, baremetal_products, _ = parse_reaction(reaction_baremetal)
            _, theo_products, _ = parse_reaction(reaction_theo)

            baremetal_reac_set = set()
            for compound in baremetal_products.items():
                baremetal_reac_set.add(compound)
            
            theo_reac_set = set()
            for compound in theo_products.items():
                theo_reac_set.add(compound)

            reac_combination = [theo_reac_set, baremetal_reac_set]
            if reac_combination not in reactions: 
                reactions.append(reac_combination)
                labelval = f" normal: {make_reaction_subscripted(reaction_theo)} \n no bare metal: {make_reaction_subscripted(reaction_baremetal)}"
                label_vals.append(labelval)
            else: 
                idx = reactions.index(reac_combination)
                labelval = label_vals[idx]
            
            print(f"labelval: {labelval}")
            print(f"energy: {energy_difference[i][1][j] - energy_difference[i][0][j]}")
            color_i = np.where(product_combinations == reac_combination)[0][0]
            bar_handle = ax.bar(mu_vals[(len(mu_vals)-1-j_reversed)], (energy_difference[i][1][j] - energy_difference[i][0][j]), alpha = 1, width = round(mu_vals[1] - mu_vals[0], 2), color = scalarMap.to_rgba(color_i), label=f"{labelval}")
            
            if (bar_handle not in plot_handles): 
                plot_handles[labelval] = bar_handle
        
    subscripted_compound = electrolyte_compounds[i].translate(SUB)
    k_BT_roomtemp = 300 * 8.6173e-5 
    extended_mu_vals = np.linspace(np.min(mu_vals)-0.25, np.max(mu_vals)+0.25, 50)
    ax.plot(extended_mu_vals, (np.ones_like(extended_mu_vals) * k_BT_roomtemp), color = "red", linestyle = "--")
    t = ax.text(np.min(mu_vals)-0.15, k_BT_roomtemp, f"kBT at 300K", color = "red")
    t.set_bbox(dict(facecolor='white', alpha=1, edgecolor='white'))
    ax.set_title(subscripted_compound)
    ax.set_xlim(np.min(mu_vals)-0.25,np.max(mu_vals)+0.25)
    ax.set_ylabel("Difference in reaction \ngrand potentials \n(eV/atom)")
    ax.set_xlabel("Potential of Li or Na (eV)")
    fig.savefig(f"../{plot_dirpath}/reaction_theo_vs_nobaremetal_{electrolyte_compounds[i]}.png",dpi=300)
    #legend = ax.legend(plot_handles.values(), plot_handles.keys(), bbox_to_anchor=(-1,-1), loc="center", fontsize = "small")
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
    # Put a legend below current axis
    legend = ax.legend(plot_handles.values(), plot_handles.keys(), loc='upper center', bbox_to_anchor=(0.5, -0.2))
    legend.get_frame().set_alpha(1)       
    figure = plt.gcf() # get current figure
    figure.set_size_inches(12, 6)
    fig.savefig(f"../{plot_dirpath}/reaction_theo_vs_nobaremetal_{electrolyte_compounds[i]}_LEGEND.png",dpi=300, bbox_inches='tight')
    ax.clear()