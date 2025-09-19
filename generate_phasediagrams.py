#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 09:59:17 2025

@author: sgreene
"""

print("loading packages")

#from mp_api.client import MPRester
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
import matplotlib.pyplot as plt
import pandas 
from collections import OrderedDict

def get_entry(comp, entries):

    for entry in entries:
        
        if (entry.reduced_formula == comp): return entry
    
    return None


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

    #print(f"react_line: {react_line}")
    ### initialize dictionaries for reactants and products ###
    reactant_dict = {}
    product_dict = {}

    ### get energy ###
    #print(f"react_line pre: {react_line}")
    react_energy = float(react_line.split()[-2])
    react_line = " ".join(react_line.split()[:-2])
    #print(f"react_line post: {react_line}")
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


def generate_GPPD_and_GPIR(elements_, relative_mu, composition_, inputfile=None, comp_benign=None, open_el=None, input_products=None):
    mpr = MPRester('GV69Do38gPOc4uk2Y4ZjFk9Wwg9p5xgw')

    '''
    loading in entries from materials project as ComputedStructureEntries
    '''
    ### elements__ = [li, La, Zr, O]
    
    ### grabs all entries from materials project that permutations of the ###
    ### provided elements and satisfy the additional criteria ###
    entries = mpr.get_entries_in_chemsys(elements_) #, additional_criteria={"thermo_types":["R2SCAN"]})
    #entries =  mpr.get_entries_in_chemsys(elements_,
    pd_entries = entries
    for entry in pd_entries: 
        if (len(entry.composition.formula.split()) == 1) and (open_el in entry.composition.formula):
            #print(f"original energy: {entry._energy}")
            entry._energy += relative_mu
            #print(f"entry.composition: {entry.composition.formula}  {entry.composition.formula.split()}")
            #print(f"final energy: {entry._energy}\n\n")


    '''
    for entry in entries: 
        print(f"entry.composition: {entry.composition}")
        print(f"entry.energy: {entry.energy}")
    '''
    
    print("vasprun init step")

    '''
    settings for loading in vasprun object
    '''
    if (inputfile == None): pass
    else: 
        filename = inputfile            # xml file name (str)
        ionic_step_skip = None          # read every "ith" frame if set to "i" (int, default = None))
        ionic_step_offset = -1          # start reading at "ith" frame if set to "i" (int, default = 0)
        parse_dos = False               # read in density of states (bool, default = True))
        parse_eigen = False             # read in eigenvalues (bool, default = True))
        parse_projected_eigen = False   # read in projected eigenvals / magnetization (bool, default = True))
        parse_potcar_file = False       # read in potcar (bool, default = True))
        occu_tol = 1e-8                 # minimum tolerance of vbm and cbm (float, default = 1e-8)
        separate_spins= False           # report vbm, cbm, and band gap for each individual spin channel (bool, default = False, must be spin-polarized calc.)
        exception_on_bad_xml = True     # throw exception on bad parsing (bool, default = False)

        vasprun_obj = Vasprun(filename, ionic_step_skip, ionic_step_offset, parse_dos, parse_eigen, parse_projected_eigen, 
                    parse_potcar_file, occu_tol, separate_spins, exception_on_bad_xml)


        '''
        converting loaded Vaspun object to ComputedStructureEntry object
        '''
        #computed entry object, grab last structure's energy from relaxation

        inc_structure = True          # read in ComputedStructureEntry instead of ComputedEntry (bool, defalut = True)
        parameters = None              # input parameters supported by Vasprun object (list, default = None, default set of params for post-processing given)
        data = None                    # output data to include supported by the Vasprun object (dict, default = None)
        entry_id = None                # entry id for the ComputedEntry (str, default = "vasprun-{current datetime}”)

        CSE_loaded =  vasprun_obj.get_computed_entry(inc_structure, parameters, data, entry_id)
        #print(f"type(CSE_loaded): {type(CSE_loaded)}")
        entries.append(CSE_loaded)


    if (input_products == None): pass
    else: 
        ionic_step_skip = None          # read every "ith" frame if set to "i" (int, default = None))
        ionic_step_offset = -1          # start reading at "ith" frame if set to "i" (int, default = 0)
        parse_dos = False               # read in density of states (bool, default = True))
        parse_eigen = False             # read in eigenvalues (bool, default = True))
        parse_projected_eigen = False   # read in projected eigenvals / magnetization (bool, default = True))
        parse_potcar_file = False       # read in potcar (bool, default = True))
        occu_tol = 1e-8                 # minimum tolerance of vbm and cbm (float, default = 1e-8)
        separate_spins= False           # report vbm, cbm, and band gap for each individual spin channel (bool, default = False, must be spin-polarized calc.)
        exception_on_bad_xml = True     # throw exception on bad parsing (bool, default = False)
        
        for product_filename in input_products:
        
            vasprun_obj = Vasprun(product_filename, ionic_step_skip, ionic_step_offset, parse_dos, parse_eigen, parse_projected_eigen, 
                        parse_potcar_file, occu_tol, separate_spins, exception_on_bad_xml)

            '''
            converting loaded Vaspun object to ComputedStructureEntry object
            '''
            #computed entry object, grab last structure's energy from relaxation

            inc_structure = True          # read in ComputedStructureEntry instead of ComputedEntry (bool, defalut = True)
            parameters = None              # input parameters supported by Vasprun object (list, default = None, default set of params for post-processing given)
            data = None                    # output data to include supported by the Vasprun object (dict, default = None)
            entry_id = None                # entry id for the ComputedEntry (str, default = "vasprun-{current datetime}”)

            CSE_loaded =  vasprun_obj.get_computed_entry(inc_structure, parameters, data, entry_id)
            #print(f"type(CSE_loaded): {type(CSE_loaded)}")
            print(f"CSE_loaded: energy: {CSE_loaded.energy}    composition: {CSE_loaded.composition}")
            entries.append(CSE_loaded)
    

    print(f"len(entries): {len(entries)}")
    print(f"entries: {entries}")
    print(f"entries[0].composition: {entries[4].reduced_formula}")

    comps_screened = []
    entries_screened = []
    desired_comps = ["Li13La6Zr3TaO24", "Zr", "Ta", "Li", "La", "LiLa2TaO6", "Li6Zr2O7", 
                    "La2O3", "Zr4O", "Li2O", "Li7La3Zr2O12"]

    for entry in entries: 
        print(entry.reduced_formula)

        entry_found = False
        if (entry.reduced_formula in desired_comps):
            print("accepted")
            if not entries_screened:    
                entries_screened.append(entry)
            else: 
                for i in range(len(entries_screened)): 
                    documented_entry = entries_screened[i]
                    if (documented_entry.reduced_formula == entry.reduced_formula):
                        entry_found = True 
                        if (entry.energy < documented_entry.energy): 
                            entries_screened.pop(i)
                            entries_screened.append(entry)

            if (entry_found == False): 
                entries_screened.append(entry)
                comps_screened.append(entry.reduced_formula)
    
    print(f"entries_screened: {entries_screened}")
    print(f"comps_screened: {comps_screened}")
    #entries = entries_screened
    # in case of GGA+U calculations
    #compat = MaterialsProject2020Compatibility()
    #entries = compat.process_entries(entries)
   
   

    pd = PhaseDiagram(pd_entries)
    
    print(f"PD Stable Entries:")
    for e in pd.stable_entries: 
        print(f"PD stable: {e.composition}")
    
    print(f"\n")
    print(f"PD Unstable Entries:")
    for e in pd.unstable_entries: 
        print(f"PD unstable: {e.composition}")


    print(f"\n")
    
    # Get the chemical potential of the pure subtance.
    mu = pd.get_transition_chempots(Element(open_el))[0]
    
    # Set the chemical potential in the elemental reservoir.
    chempots = {open_el: relative_mu + mu}
    

    # Build the grand potential phase diagram
    gpd = GrandPotentialPhaseDiagram(entries=entries, chempots=chempots) #, elements=elements_)
    
    print(f"\n")
    print(f"\nGPD Stable Entries:")
    for e in gpd.stable_entries: 
        print(f"GPPD stable: {e.composition}")


    print(f"\n")
    print(f"GPD Unstable Entries:")
    for e in gpd.unstable_entries: 
        print(f"GPPD unstable: {e.composition}")

    comp_SE = Composition(composition_)

    GPIR_reacs = GrandPotentialInterfacialReactivity(c1=comp_SE, c2=comp_benign, grand_pd=gpd, pd_non_grand=pd, norm=True, include_no_mixing_energy=True, use_hull_energy=False)
    

    critical_rxns = [
        OrderedDict([
            ("Atomic fraction", round(ratio, 4)),
            ("Reaction equation", rxn),
            ("E$_{rxt}$ per mol equation (kJ/mol)", round(rxn_energy, 4)),
            ("E$_{rxt}$ per reactant atom (eV/atom)", round(reactivity, 4)),
        ])

        for _, ratio, reactivity, rxn, rxn_energy in GPIR_reacs.get_kinks()]
    interface_reaction_table = pandas.DataFrame(critical_rxns)

    with pandas.option_context('display.max_rows', None, 'display.max_columns', None, 'display.max_colwidth', None):  # more options can be specified also
        print(interface_reaction_table)

    
    allatoms_str = '_'.join(elements_)
    
    comp_SE = Composition(composition_)

    IR_reacs = InterfacialReactivity(c1=comp_SE, c2=Composition(open_el), pd=pd,  norm=True, include_no_mixing_energy=True, use_hull_energy=False)

    print(f"\n\n")

    print(f"Solid Electrolyte: {comp_SE}")

    print(f"Phase Diagram:")
    critical_rxns = [
        OrderedDict([
            ("Atomic fraction", round(ratio, 4)),
            ("Reaction equation", rxn),
            ("E$_{rxt}$ per mol equation (kJ/mol)", round(rxn_energy, 4)),
            ("E$_{rxt}$ per reactant atom (eV/atom)", round(reactivity, 4)),

        ])
        for _, ratio, reactivity, rxn, rxn_energy in IR_reacs.get_kinks()]
    interface_reaction_table = pandas.DataFrame(critical_rxns)

    with pandas.option_context('display.max_rows', None, 'display.max_columns', None, 'display.max_colwidth', None):  # more options can be specified also
        print(interface_reaction_table)

    GPIR_reacs = GrandPotentialInterfacialReactivity(c1=comp_SE, c2=comp_benign, grand_pd=gpd, pd_non_grand=pd, norm=True, include_no_mixing_energy=True, use_hull_energy=False)
    
    print(f"GPIR_reacs.gpd: {type(GPIR_reacs.pd)}")
    
    scaled_reac_es = []
    reacs = []

    for label, reac in IR_reacs.labels.items():
        #print(f"reac: {reac}")
        reac_split = reac.split()
        split_idx = 0
        while (reac_split[split_idx] != "eV/atom"): split_idx += 1
        reac_split = reac_split[(split_idx+1):]
        energy = reac_split[1]
        reac_split = reac_split[2:]
        reac_split[-1] = reac_split[-1].strip("\n")
        reac_line = " ".join(reac_split)
        reac_line = reac_line + "  " + str(energy) + " eV/atom \n"
        
        reactant_dict, product_dict, reaction_energy = parse_reaction(reac_line)

        Li_coeff = 0
        total_coeff = 0

        for whole_compound, whole_coeff in reactant_dict.items():
            
            compound_split = seperate_string_number(whole_compound)
            coefficients = np.asarray([int(num) for num in compound_split if num.isdigit()])
            elements = np.asarray([str(elem) for elem in compound_split if not elem.isdigit()])

            for i in range(len(coefficients)): 
                if (elements[i] == open_el): Li_coeff += coefficients[i] * whole_coeff
                total_coeff += coefficients[i] * whole_coeff
            

        if ( (total_coeff - Li_coeff) != 0):
            scaled_reac_energy = reaction_energy *  total_coeff / (total_coeff - Li_coeff) 
            scaled_reac_es.append(scaled_reac_energy)
            reacs.append(reac_line)

    scaled_reac_es = np.asarray(scaled_reac_es)
    min_energy = np.min(scaled_reac_es)
    min_idx = np.where(scaled_reac_es == min_energy)[0][0]
    min_reac = reacs[min_idx]

    critical_rxns = [
        OrderedDict([
            ("Atomic fraction", round(ratio, 4)),
            ("Reaction equation", rxn),
            ("E$_{rxt}$ per mol equation (kJ/mol)", round(rxn_energy, 4)),
            ("E$_{rxt}$ per reactant atom (eV/atom)", round(reactivity, 4)),
        ])

        for _, ratio, reactivity, rxn, rxn_energy in GPIR_reacs.get_kinks()]
    interface_reaction_table = pandas.DataFrame(critical_rxns)

    with pandas.option_context('display.max_rows', None, 'display.max_columns', None, 'display.max_colwidth', None):  # more options can be specified also
        print(interface_reaction_table)
    
    '''
    print(f"Grand Potential Phase Diagram no-mix:")
    critical_rxns_nomix = [
        OrderedDict([
            ("Atomic fraction", round(ratio, 4)),
            ("Reaction equation", rxn),
            ("E$_{rxt}$ per mol equation (kJ/mol)", round(rxn_energy, 4)),
            ("E$_{rxt}$ per reactant atom (eV/atom)", round(reactivity, 4)),
        ])

        for _, ratio, reactivity, rxn, rxn_energy in GPIR_reacs_nomix.get_kinks()]
    interface_reaction_table = pandas.DataFrame(critical_rxns_nomix)

    with pandas.option_context('display.max_rows', None, 'display.max_columns', None, 'display.max_colwidth', None):  # more options can be specified also
        print(interface_reaction_table)
    '''
    
    allatoms_str = '_'.join(elements)
    react_output = "all_reactions/" + allatoms_str + f"_{round(relative_mu, 3)}volts_reactions.txt"
    
    react_file = open(react_output, "w")

    react_file.write(f"{composition_} decomposition reactions\n--------\n")


    for key, value in GPIR_reacs.labels.items():
        split_line = value.split()
        # print(value)
        parse_i = 0
        while ("eV/atom" not in split_line[parse_i]): parse_i += 1
        react_energy = split_line[(parse_i + 2)]
        split_value = split_line[(parse_i + 3):]
        join_value = " ".join(split_value) + f"  {react_energy} eV/atom"
        react_file.write(f"{join_value}\n")


    react_file.close()

    print("grand potential phase diagram has been made!")
    

def generate_PD(elements_, inputfile=None, input_products=None):
    mpr = MPRester('JxCBb3dqR7jtlimd8RA9DokS7uX8Ru7m')

    '''
    loading in entries from materials project as ComputedStructureEntries
    '''
    ### elements__ = [li, La, Zr, O]
    
    ### grabs all entries from materials project that permutations of the ###
    ### provided elements and satisfy the additional criteria ###
    entries = mpr.get_entries_in_chemsys(elements_) #, additional_criteria={"thermo_types":["R2SCAN"]})
    #entries =  mpr.get_entries_in_chemsys(elements_,
    pd_entries = entries

    
    print("vasprun init step")

    '''
    settings for loading in vasprun object
    '''
    if (inputfile == None): pass
    else: 
        filename = inputfile            # xml file name (str)
        ionic_step_skip = None          # read every "ith" frame if set to "i" (int, default = None))
        ionic_step_offset = -1          # start reading at "ith" frame if set to "i" (int, default = 0)
        parse_dos = False               # read in density of states (bool, default = True))
        parse_eigen = False             # read in eigenvalues (bool, default = True))
        parse_projected_eigen = False   # read in projected eigenvals / magnetization (bool, default = True))
        parse_potcar_file = False       # read in potcar (bool, default = True))
        occu_tol = 1e-8                 # minimum tolerance of vbm and cbm (float, default = 1e-8)
        separate_spins= False           # report vbm, cbm, and band gap for each individual spin channel (bool, default = False, must be spin-polarized calc.)
        exception_on_bad_xml = True     # throw exception on bad parsing (bool, default = False)

        vasprun_obj = Vasprun(filename, ionic_step_skip, ionic_step_offset, parse_dos, parse_eigen, parse_projected_eigen, 
                    parse_potcar_file, occu_tol, separate_spins, exception_on_bad_xml)


        '''
        converting loaded Vaspun object to ComputedStructureEntry object
        '''
        #computed entry object, grab last structure's energy from relaxation

        inc_structure = True          # read in ComputedStructureEntry instead of ComputedEntry (bool, defalut = True)
        parameters = None              # input parameters supported by Vasprun object (list, default = None, default set of params for post-processing given)
        data = None                    # output data to include supported by the Vasprun object (dict, default = None)
        entry_id = None                # entry id for the ComputedEntry (str, default = "vasprun-{current datetime}”)

        CSE_loaded =  vasprun_obj.get_computed_entry(inc_structure, parameters, data, entry_id)
        #print(f"type(CSE_loaded): {type(CSE_loaded)}")
        entries.append(CSE_loaded)

    if (input_products == None): pass
    else: 
        ionic_step_skip = None          # read every "ith" frame if set to "i" (int, default = None))
        ionic_step_offset = -1          # start reading at "ith" frame if set to "i" (int, default = 0)
        parse_dos = False               # read in density of states (bool, default = True))
        parse_eigen = False             # read in eigenvalues (bool, default = True))
        parse_projected_eigen = False   # read in projected eigenvals / magnetization (bool, default = True))
        parse_potcar_file = False       # read in potcar (bool, default = True))
        occu_tol = 1e-8                 # minimum tolerance of vbm and cbm (float, default = 1e-8)
        separate_spins= False           # report vbm, cbm, and band gap for each individual spin channel (bool, default = False, must be spin-polarized calc.)
        exception_on_bad_xml = True     # throw exception on bad parsing (bool, default = False)
        
        for product_filename in input_products:
        
            vasprun_obj = Vasprun(product_filename, ionic_step_skip, ionic_step_offset, parse_dos, parse_eigen, parse_projected_eigen, 
                        parse_potcar_file, occu_tol, separate_spins, exception_on_bad_xml)

            '''
            converting loaded Vaspun object to ComputedStructureEntry object
            '''
            #computed entry object, grab last structure's energy from relaxation

            inc_structure = True          # read in ComputedStructureEntry instead of ComputedEntry (bool, defalut = True)
            parameters = None              # input parameters supported by Vasprun object (list, default = None, default set of params for post-processing given)
            data = None                    # output data to include supported by the Vasprun object (dict, default = None)
            entry_id = None                # entry id for the ComputedEntry (str, default = "vasprun-{current datetime}”)

            CSE_loaded =  vasprun_obj.get_computed_entry(inc_structure, parameters, data, entry_id)
            #print(f"type(CSE_loaded): {type(CSE_loaded)}")
            print(f"CSE_loaded: energy: {CSE_loaded.energy}    composition: {CSE_loaded.composition}")
            entries.append(CSE_loaded)
    

    # in case of GGA+U calculations
    #compat = MaterialsProject2020Compatibility()
    #entries = compat.process_entries(entries)   

    pd = PhaseDiagram(pd_entries)
    
    print(f"PD Stable Entries:")
    for e in pd.stable_entries: 
        print(f"PD stable: {e.composition}")
    
    print(f"\n")
    print(f"PD Unstable Entries:")
    for e in pd.unstable_entries: 
        print(f"PD unstable: {e.composition}")


    print(f"\n")

    # Plot phase diagram
    PDPlotter(pd).show()
    

mu_vals = np.asarray([0]) 

for mu_in in mu_vals:
    print(f"mu_in {mu_in} volts")
    input_filename = "vasprun_LLTZO.xml" 
    inputproducts = ["vasprun_LiLa2TaO6.xml"]
    elements = ['Li', 'La', 'O', 'Ta'] # LLTO, can only use 4 elements if you want to visualize it
    generate_PD(elements, input_products=inputproducts) #, inputfile = input_filename) # 
    
    
    print(f"mu_in {mu_in} volts")
    input_filename = "vasprun_LLTZO.xml" 
    inputproducts = ["vasprun_LiLa2TaO6.xml"]
    elements = ['Li', 'La', 'O', 'Ta'] # LLZTO, no visualization
    generate_GPPD_and_GPIR(elements, mu_in, "Li13O24La6Zr3Ta1", inputfile = input_filename , comp_benign=Composition("Li2O"), open_el="Li", input_products=inputproducts)
    