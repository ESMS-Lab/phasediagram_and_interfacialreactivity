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


def generate_GPPD_and_GPIR(elements_, relative_mu, composition_, write_dir = None, inputfile=None, comp_benign=None, open_el=None, input_products=None, only_exp_observed=False, unwanted_elems=[]):
    mpr = MPRester('JxCBb3dqR7jtlimd8RA9DokS7uX8Ru7m')

    ### grabs all entries from materials project that permutations of the ###
    ### provided elements and satisfy the additional criteria ###
    
    if (only_exp_observed):
        # First obtain all possible subsets of the chemical system
        constituent_chemsys = set()
        for i in range(len(elements_)):
            constituent_chemsys.update({"-".join(combo) for combo in combinations(elements_,1+i)})

        # This is the default thermo type used in MP, you could also use GGA_GGA+U or R2SCAN.
        thermo_types = ["GGA_GGA+U_R2SCAN"]

        # with MPRester("your_api_key") as mpr:
        # Retrieve all materials in those chemical systems which are experimentally observed, `theoretical=False`
        expt_obs_mats = mpr.materials.summary.search(chemsys = list(constituent_chemsys), theoretical=False, fields=["material_id","composition"])
        
        # Obtain thermodynamic data
        thermo_docs = mpr.materials.thermo.search(material_ids=[doc.material_id for doc in expt_obs_mats],thermo_types=thermo_types)

        # Concatenate entries in the thermo docs
        exp_entries = []
        for doc in thermo_docs:
            exp_entries.extend(doc.entries.values())

    theo_exp_entries = mpr.get_entries_in_chemsys(elements_) #, additional_criteria={"thermo_types":["R2SCAN"]})

    entries = []
    bool_array = np.ones_like(theo_exp_entries).astype(bool)
    
    if (only_exp_observed):
        for i in range(len(exp_entries)):
            exp_entry = exp_entries[i]
            #print(exp_entry.formula)
            entry_found = False
            #print(f"exp_entry test: {exp_entry.composition}")
            for j in range(len(theo_exp_entries)): 
                theo_exp_entry = theo_exp_entries[j] 
                
                if (exp_entry.formula == theo_exp_entry.formula):
                    #print(f"FOUND ENTRY: {theo_exp_entry.composition}")
                    bool_array[j] = False
                    if ((len(unwanted_elems) != 0) and (len(theo_exp_entry.elements) == 1)):
                        #print(f"theo_exp_entry.elements first if statement: {theo_exp_entry.elements}")
                        if (theo_exp_entry.elements[0].symbol in unwanted_elements):
                            print(theo_exp_entry.composition)
                            print(f"old energy: {theo_exp_entry._energy}")
                            theo_exp_entry._energy += 1000
                            print(f"new energy: {theo_exp_entry._energy}")

                    entries.append(theo_exp_entry)
                    entry_found = True
                    
            if not entry_found:
                print(f"ENTRY NOT FOUND: {exp_entry.composition} energy: {exp_entry.energy}")
                entries.append(exp_entry)
            
                    
    else: entries = theo_exp_entries
    

    '''
    settings for loading in dft-calculated structures not present in Materials Project
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
            #print(f"CSE_loaded: energy: {CSE_loaded.energy}    composition: {CSE_loaded.composition}")
            entries.append(CSE_loaded)
    

    pd_entries = entries
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
    print(f"GPPD Stable Entries:")
    for e in gpd.stable_entries: 
        print(f"GPPD stable: {e.composition}")
    
    print(f"\n")
    print(f"GPPD Unstable Entries:")
    
    for e in gpd.unstable_entries: 
        print(f"GPPD unstable: {e.composition}")
    

    comp_SE = Composition(composition_)
    GPIR_reacs = GrandPotentialInterfacialReactivity(c1=comp_SE, c2=comp_benign, grand_pd=gpd, pd_non_grand=pd, norm=True, include_no_mixing_energy=True, use_hull_energy=False)
    print(f"comp_SE: {comp_SE}")
    print(f"comp_benign: {comp_benign}")
    
    print(f"GPIR_reacs.comp1: {GPIR_reacs.comp1}")
    print(f"GPIR_reacs.comp2: {GPIR_reacs.comp2}")
    critical_rxns = [
        OrderedDict([
            ("Atomic fraction", round(ratio, 4)),
            ("Reaction equation", rxn),
            ("E$_{rxt}$ per mol equation (kJ/mol)", round(rxn_energy, 4)),
            ("E$_{rxt}$ per reactant atom (eV/atom)", round(reactivity, 4)),
        ])

        for _, ratio, reactivity, rxn, rxn_energy in GPIR_reacs.get_kinks()]
    GPIR_interface_reaction_table = pandas.DataFrame(critical_rxns)
    with pandas.option_context('display.max_rows', None, 'display.max_columns', None, 'display.max_colwidth', None):  # more options can be specified also
        print(f"GPIR_interface_reaction_table: {GPIR_interface_reaction_table}")
    
    if (write_dir != None):
        allatoms_str = '_'.join(elements)
        react_output = allatoms_str + f"_{round(relative_mu, 3)}volts_reactions.txt"

        pwd = os.getcwd()
        react_dir = os.path.join(pwd, "no_baremetal_reactions")
        #react_dir = os.path.join(pwd, "all_reactions")

        if os.path.isdir(react_dir): pass
        else: os.mkdir(react_dir)
        react_output = os.path.join(react_dir, react_output)
        print(f"react_dir: {react_dir}")
        print(f"react_output: {react_output}")
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

        '''
        writing stable compounds to output files 
        '''
        allatoms_str = '_'.join(elements)
        
        decomp_output = allatoms_str + f"_{round(relative_mu, 3)}volts_decomp_output.txt"
        if (write_dir != None):
            pwd = os.getcwd()
            decomp_dir = os.path.join(pwd, write_dir)
            if os.path.isdir(decomp_dir): pass
            else: os.mkdir(decomp_dir)
            
            decomp_output = os.path.join(decomp_dir, decomp_output)
        
        decomp_dir = os.path.join(pwd, write_dir)    
        stable_output = allatoms_str + "_stable_decomp_output.txt"
        stable_output = os.path.join(decomp_dir, stable_output)
        stableonly_file = open(stable_output, "a")
        stableonly_file.write(f"Stable Entries (formula) {round(relative_mu, 3)}volts\n--------\n")
        
        for e in gpd.stable_entries:
            stableonly_file.write(f"{' '.join(str(e).split()[4:])}\n")

        stableonly_file.write(f"\n")
        stableonly_file.close()   


def generate_PD_and_IR(elements_, reactant1, reactant2, input_reactants=None, input_products=None, only_exp_observed=False, unwanted_elems=[]):
    mpr = MPRester('JxCBb3dqR7jtlimd8RA9DokS7uX8Ru7m')

    '''
    loading in entries from materials project as ComputedStructureEntries
    '''
    mpr = MPRester('JxCBb3dqR7jtlimd8RA9DokS7uX8Ru7m')

    ### grabs all entries from materials project that permutations of the ###
    ### provided elements and satisfy the additional criteria ###
    
    if (only_exp_observed):
        # First obtain all possible subsets of the chemical system
        constituent_chemsys = set()
        for i in range(len(elements_)):
            constituent_chemsys.update({"-".join(combo) for combo in combinations(elements_,1+i)})

        #print(f"constituent_chemsys: {constituent_chemsys}")
        # This is the default thermo type used in MP, you could also use GGA_GGA+U or R2SCAN.
        thermo_types = ["GGA_GGA+U_R2SCAN"]

        # with MPRester("your_api_key") as mpr:
        # Retrieve all materials in those chemical systems which are experimentally observed, `theoretical=False`
        expt_obs_mats = mpr.materials.summary.search(chemsys = list(constituent_chemsys), theoretical=False, fields=["material_id","composition"])
        
        #print(f"expt_obs_mats: {expt_obs_mats}")

        # Obtain thermodynamic data
        thermo_docs = mpr.materials.thermo.search(material_ids=[doc.material_id for doc in expt_obs_mats],thermo_types=thermo_types)

        # Concatenate entries in the thermo docs
        exp_entries = []
        for doc in thermo_docs:
            exp_entries.extend(doc.entries.values())

    theo_exp_entries = mpr.get_entries_in_chemsys(elements_) #, additional_criteria={"thermo_types":["R2SCAN"]})

    entries = []
    bool_array = np.ones_like(theo_exp_entries).astype(bool)
    
    if (only_exp_observed):
        for i in range(len(exp_entries)):
            exp_entry = exp_entries[i]
            entry_found = False
            for j in range(len(theo_exp_entries)): 
                theo_exp_entry = theo_exp_entries[j] 
                
                if (exp_entry.formula == theo_exp_entry.formula):
                    bool_array[j] = False
                    if ((len(unwanted_elems) != 0) and (len(theo_exp_entry.elements) == 1)):
                        #print(f"theo_exp_entry.elements first if statement: {theo_exp_entry.elements}")
                        if (theo_exp_entry.elements[0].symbol in unwanted_elements):
                            print(theo_exp_entry.composition)
                            print(f"old energy: {theo_exp_entry._energy}")
                            theo_exp_entry._energy += 1000
                            print(f"new energy: {theo_exp_entry._energy}")

                    entries.append(theo_exp_entry)
                    entry_found = True
                    
            if not entry_found:
                print(f"ENTRY NOT FOUND: {exp_entry.composition} energy: {exp_entry.energy}")
                entries.append(exp_entry)
            
                    
    else: entries = theo_exp_entries
    
    theo_exp_entries = np.asarray(theo_exp_entries)
    masked_theo_exp = theo_exp_entries[bool_array]
    
    print(f"screened only theo_exp_entries: {set([entry.reduced_formula for entry in masked_theo_exp])}")
        
    print("vasprun init step")

    '''
    settings for loading in dft-calculated structures not present in Materials Project
    '''
    if (input_reactants == None): pass
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
            #print(f"CSE_loaded: energy: {CSE_loaded.energy}    composition: {CSE_loaded.composition}")
            entries.append(CSE_loaded)
    

    exp_entries_set = np.sort(list((set([tuple([entry.composition.reduced_formula, entry.energy_per_atom]) for entry in entries]))))
    theo_exp_entries_set = np.sort(list((set([tuple([entry.composition.reduced_formula, entry.energy_per_atom]) for entry in theo_exp_entries]))))

    print(f"exp_entries_set: {exp_entries_set}")
    print(f"theo_exp_entries_set: {theo_exp_entries_set}")
    
    print(f"only in exp_entries_set:")
    for item in exp_entries_set:
        if item not in theo_exp_entries_set:
            print(item)

    print(f"only in theo_exp_entries:")
    for item in theo_exp_entries_set:
        if item not in exp_entries_set:
            print(item)

    pd_entries = entries
    
    print("vasprun init step")

    '''
    settings for loading in vasprun object
    '''
    if (input_reactants == None): pass
    else: 
        filename = input_reactants            # xml file name (str)
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
        filename = input_products
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

    if (len(Composition(reactant1).elements) > len(elements_) ) or (len(Composition(reactant2).elements) > len(elements_) ):
        print(f"Cannot use InterfacialReactivity modlue, either {reactant1} or {reactant2} has elements not in {elements}")

    else: 
        ### generating interfacial reactions ###
        IR_reacs = InterfacialReactivity(c1=Composition(reactant1), c2=Composition(reactant2), pd=pd,  norm=True, include_no_mixing_energy=True, use_hull_energy=False)

        print(f"\n\n")
        print(f"Reactant 1: {Composition(reactant1)}")
        print(f"Reactant 2: {Composition(reactant2)}")

        print(f"Phase Diagram:")
        critical_rxns = [
            OrderedDict([
                ("Atomic fraction", round(ratio, 4)),
                ("Reaction equation", rxn),
                ("E$_{rxt}$ per mol equation (kJ/mol)", round(rxn_energy, 4)),
                ("E$_{rxt}$ per reactant atom (eV/atom)", round(reactivity, 4)),

            ])
            for _, ratio, reactivity, rxn, rxn_energy in IR_reacs.get_kinks()]
        IR_interface_reaction_table = pandas.DataFrame(critical_rxns)
        with pandas.option_context('display.max_rows', None, 'display.max_columns', None, 'display.max_colwidth', None):  # more options can be specified also
            print(f"IR_interface_reaction_table: {IR_interface_reaction_table}")
    

mu_vals = np.asarray([0]) 

for mu_in in mu_vals:
    print(f"mu_in {mu_in} volts")
    input_filename = "vasprun_LLTZO.xml" 
    inputproducts = ["vasprun_LiLa2TaO6.xml"]
    inputreactants = [] 
    only_exp_observed_val = False
    elements = ['Li', 'La', 'O', 'Ta', 'Zr'] # LLTO, can only use 4 elements if you want to visualize it
    reactant_1 = "Li13O24La6Zr3Ta1"
    reactant_2 = "Li2O"
    generate_PD_and_IR(elements, reactant_1, reactant_2, input_products=inputproducts, only_exp_observed=only_exp_observed_val) #, inputfile = input_filename) # 
    
    
    print(f"mu_in {mu_in} volts")
    input_filename = "vasprun_LLTZO.xml" 
    inputproducts = ["vasprun_LiLa2TaO6.xml"]
    elements = ['Li', 'La', 'O', 'Ta', 'Zr'] # LLZTO, no visualization
    generate_GPPD_and_GPIR(elements, mu_in, "Li13O24La6Zr3Ta1", inputfile = input_filename , comp_benign=Composition("Li2O"), open_el="Li", input_products=inputproducts)
    