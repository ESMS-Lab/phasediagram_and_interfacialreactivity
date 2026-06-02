from mp_api.client import MPRester
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
import pandas 
from collections import OrderedDict
from itertools import combinations


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


def get_min_energy_entry(mat_comp_str, entries):

    with MPRester(api_key="JxCBb3dqR7jtlimd8RA9DokS7uX8Ru7m") as mpr:
        energy_min = sys.maxsize
        entry_min = None    
        mat_comp = Composition(mat_comp_str)
        element_of_interest = mat_comp.elements[0]
        mat_comp_dict = mat_comp.get_el_amt_dict()
        mat_comp_amt = mat_comp_dict[element_of_interest.symbol]

        ### iterating through different entries in MP ###
        for entry in entries:
            #print(f"entry.reduced_formula: {entry.reduced_formula}")
            if (entry.reduced_formula == mat_comp_str):
                value_comp_dict = entry.composition.get_el_amt_dict()
                value_amt = value_comp_dict[element_of_interest.symbol]
                scale_factor = mat_comp_amt / value_amt
                energy = entry.energy * scale_factor
                #print(f"accepted formula: {entry.reduced_formula}")

                if (energy < energy_min): 
                    energy_min = energy
                    entry_min = entry

                if (entry_min == None): 
                    ValueError(f"MP entry objects being incorrectly accessed or requested material does not exist")
                    exit()
        
    return energy_min, entry_min


def get_num_atoms(compound):
    compound_split = seperate_string_number(compound)
    coefficients = np.asarray([int(num) for num in compound_split if num.isdigit()])
    elements = np.asarray([str(elem) for elem in compound_split if not elem.isdigit()])
    
    num_atoms = 0
    num_nonopen_atoms = 0
    
    #print(f"elements: {elements}  coefficients: {coefficients}")
    for i in range(len(elements)):
        if not (elements[i] in ["Li", "Na"]): num_nonopen_atoms += coefficients[i]
        num_atoms += coefficients[i]

    #print(f"num_nonopen_atoms: {num_nonopen_atoms}")
    return num_atoms, num_nonopen_atoms


def is_number_or_dec(s):
    """ Returns True if string is a number. """
    return s.replace('.','',1).isdigit()


def read_select_theoretical_compounds_file(filename, elements_in_react):
    file = open(filename, 'r')
    lines = [line for line in file]
    selected_compounds = []
    
    # print(f"elements_in_react: {elements_in_react}")
    for line in lines:
        split_line = line.split()
        
        if (split_line[0] == "X") or (split_line[0] == "?"): 
            pass
        else: 
            compound = Composition(split_line[0].strip(","))
            elems_dict = compound.get_el_amt_dict()
            # print(f"compoud: {compound}")

            compound_in_subspace = True
            for elem, amount in elems_dict.items():
                # print(f"elem: {elem}")
                if elem not in elements_in_react: 
                    # print(f"element not in subspace: {elem}")
                    compound_in_subspace = False

            if compound_in_subspace: 
                selected_compounds.append(compound)
                # print(f"accepted custom compound: {compound}")

    selected_compounds = set(selected_compounds)

    return selected_compounds


def identify_lowest_e_entry(composition, entries):
    entries_set = np.sort(list((set([tuple([entry.composition.reduced_formula, entry.energy_per_atom]) for entry in entries]))))

    located_compoud = None

    for entry in entries_set:     
        # print(f"composition.reduced_formula: {composition.reduced_formula}  entry[0]: {entry[1]}")    
        if (composition.reduced_formula == entry[1]):
            if (located_compoud is None): located_compoud = entry
            elif (float(located_compoud[0]) > float(entry[0])): located_compoud = entry
            
    if (located_compoud is None):  
        ValueError(f"Can't find the desired compound in the set of entries")
        print(f"ERROR: Can't find the desired compound in the set of entries")
        exit()

    return located_compoud


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


def generate_GPPD_and_GPIR(elements_, relative_mu, composition_, write_dir = None, inputfile=None, comp_benign=None, open_el=None, input_products=None, only_exp_observed=False, _unwanted_elems=[], custom_theo_entries_file = None):
    mpr = MPRester('JxCBb3dqR7jtlimd8RA9DokS7uX8Ru7m')

    ### grabs all entries from materials project that permutations of the ###
    ### provided elements and satisfy the additional criteria ###

    if (only_exp_observed) or (custom_theo_entries_file != None):
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
            entry_found = False

            for j in range(len(theo_exp_entries)): 
                theo_exp_entry = theo_exp_entries[j] 
                
                if (exp_entry.reduced_formula == theo_exp_entry.reduced_formula):
                    
                    bool_array[j] = False
                    if ((len(_unwanted_elems) != 0) and (len(theo_exp_entry.elements) == 1)):
                        
                        if (theo_exp_entry.elements[0].symbol in unwanted_elements):
                            #print(theo_exp_entry.composition)
                            #print(f"old energy: {theo_exp_entry._energy}")
                            theo_exp_entry._energy += 1000
                            #print(f"new energy: {theo_exp_entry._energy}")

                    entries.append(theo_exp_entry)
                    entry_found = True
                    
            if not entry_found:
                print(f"ENTRY NOT FOUND: {exp_entry.composition} energy: {exp_entry.energy}")
                
            
    elif (custom_theo_entries_file != None): 
        #print(f"in reading file")
        for i in range(len(exp_entries)):
            exp_entry = exp_entries[i]
            entry_found = False

            for j in range(len(theo_exp_entries)): 
                theo_exp_entry = theo_exp_entries[j] 
                
                if (exp_entry.reduced_formula == theo_exp_entry.reduced_formula):
                    
                    if ((len(_unwanted_elems) != 0) and (len(theo_exp_entry.elements) == 1)):
                        
                        if (theo_exp_entry.elements[0].symbol in unwanted_elements):
                            #print(theo_exp_entry.composition)
                            #print(f"old energy: {theo_exp_entry._energy}")
                            theo_exp_entry._energy += 1000
                            #print(f"new energy: {theo_exp_entry._energy}")

                    entries.append(theo_exp_entry)
                    entry_found = True
                    
            if not entry_found:
                print(f"ENTRY NOT FOUND: {exp_entry.composition} energy: {exp_entry.energy}")

        custom_theo_entries = read_select_theoretical_compounds_file(custom_theo_entries_file, elements_)

        for entry in theo_exp_entries:
            if (entry.composition in custom_theo_entries):
                #print(f"accepted custom entry: {entry}")
                entries.append(entry)
            else: pass
                #print(f"rejected custom entry: {entry} ")

    else: entries = theo_exp_entries
    
    compat = MaterialsProject2020Compatibility()
    compat.process_entries(entries)

    '''
    exp_entries_set = np.sort(list((set([tuple([entry.composition.reduced_formula, entry.energy_per_atom]) for entry in entries]))))
    theo_exp_entries_set = np.sort(list((set([tuple([entry.composition.reduced_formula, entry.energy_per_atom]) for entry in theo_exp_entries]))))

    print(f"exp_entries_set: {exp_entries_set}")
    print(f"theo_exp_entries_set: {theo_exp_entries_set}")
    
    exp_entries_set_screened = set()
    exp_theo_entries_set_screened = set()

    print(f"exp_entries_set_screened: {exp_entries_set_screened}")
    for entry in exp_entries_set:
        compound_found = False
        
        for screened_entry in exp_entries_set_screened:
            
            if (screened_entry[1] == entry[1]):
                compound_found = True
                if (float(screened_entry[0]) > float(entry[0])):
                    #print(f"accepted entry: {entry}")
                    exp_entries_set_screened.remove(screened_entry)
                    exp_entries_set_screened.add(tuple(entry))
            
        if not compound_found:  
                
            exp_entries_set_screened.add(tuple(entry))
    
    print(f"theo_exp_entries_set: {theo_exp_entries_set}")
    for entry in theo_exp_entries_set:
        compound_found = False
        
        for screened_entry in exp_theo_entries_set_screened:
            
            if (screened_entry[1] == entry[1]):
                
                compound_found = True

                if (float(screened_entry[0]) > float(entry[0])):
                    
                    exp_theo_entries_set_screened.remove(screened_entry)
                    exp_theo_entries_set_screened.add(tuple(entry))

        if not compound_found: 
            
            exp_theo_entries_set_screened.add(tuple(entry))
        
    
    theo_exp_entries_screened = []
    print(f"only in theo_exp_entries:")
    for item in theo_exp_entries_set:
        if item not in exp_entries_set:
            print(item)
            theo_exp_entries_screened.append(item)
    
    '''
    
    write_dir = "revised_elems_nobaremetals"
    decomp_dir = os.path.join(os.getcwd(), write_dir)  
    allatoms_str = '_'.join(elements)
    react_output = allatoms_str 
    stable_output = react_output + "_theo_only_output.txt"
    stable_output = os.path.join(decomp_dir, stable_output)

    if os.path.isdir(decomp_dir): pass
    else: os.mkdir(decomp_dir)
    
    '''
    stableonly_file = open(stable_output, "a")
    stableonly_file.write(f"Stable Entries (formula) {round(relative_mu, 3)}volts\n--------\n")
    
    for item in theo_exp_entries_screened:
        stableonly_file.write(f"{item}\n")

    stableonly_file.write(f"\n")
    stableonly_file.close()   
    '''

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

        CSE_loaded_reactant =  vasprun_obj.get_computed_entry(inc_structure, parameters, data, entry_id)
        #print(f"type(CSE_loaded): {type(CSE_loaded)}")
        entries.append(CSE_loaded_reactant)

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
            
            ### altering energies of gases loaded in ###
            if ("O2_" in product_filename):
                if ("O" in elements_):
                    print(f"adding in O compounds")
                    CSE_loaded._energy = -9.57
                    entries.append(CSE_loaded)
                else: pass

            elif (f"Cl2_" in product_filename):
                if ("Cl" in elements_):
                    print(f"adding in Cl compounds")
                    CSE_loaded._energy = -4.02
                    entries.append(CSE_loaded)
                else: pass
                
            elif (f"H2_" in product_filename):
                if ("H" in elements_):
                    print(f"adding in H compounds")
                    CSE_loaded._energy = -6.99
                    entries.append(CSE_loaded)
                else: pass
                
            ### adding in all other energies ### 
            #else: 
            entries.append(CSE_loaded)
    
    SE_entry = None
    if (inputfile != None):
        SE_entry = CSE_loaded_reactant
    else:
        for entry in theo_exp_entries:
            if (entry.reduced_formula == Composition(composition_).reduced_formula):
                
                if (SE_entry is None): 
                    SE_entry = entry
                    # print(f"entry.composition: {entry.composition}")
                    entries.append(entry)
                elif (SE_entry.energy > entry.energy): 
                    SE_entry = entry
                    # print(f"entry.composition: {entry.composition}")
                    entries.append(entry)


    pd_entries = entries
    pd = PhaseDiagram(pd_entries)


    print(f"PD Stable Entries:")
    for e in pd.stable_entries: 
        print(f"PD stable: {e.composition}   energy: {e.energy}")
    
    print(f"\n")
    print(f"PD Unstable Entries:")
    for e in pd.unstable_entries: 
        print(f"PD unstable: {e.composition}  energy: {e.energy}")

    print(f"\n")
    
    # Get the chemical potential of the pure subtance.
    mu = pd.get_transition_chempots(Element(open_el))[0]
    
    # Set the chemical potential in the elemental reservoir.
    chempots = {open_el: relative_mu + mu}

    # Build the grand potential phase diagram
    gpd = GrandPotentialPhaseDiagram(entries=entries, chempots=chempots) #, elements=elements_)
    
    SE_E_above_hull = pd.get_e_above_hull(SE_entry)
    
    '''
    print(f"SE_E_above_hull: {SE_E_above_hull}")
    print(f"\n")

    print(f"GPPD Stable Entries:")
    for e in gpd.stable_entries: 
        print(f"GPPD stable: {e.composition}  energy: {e.energy}")
    
    print(f"\n")
    print(f"GPPD Unstable Entries:")
    
    for e in gpd.unstable_entries: 
        print(f"GPPD unstable: {e.composition}  energy: {e.energy}")
    '''

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
    GPIR_interface_reaction_table = pandas.DataFrame(critical_rxns)
    with pandas.option_context('display.max_rows', None, 'display.max_columns', None, 'display.max_colwidth', None):  # more options can be specified also
        print(f"GPIR_interface_reaction_table: {GPIR_interface_reaction_table}")
    
    if (write_dir != None):
        allatoms_str = '_'.join(elements)
        react_output = allatoms_str + f"_{round(relative_mu, 3)}volts_reactions.txt"

        pwd = os.getcwd()
        #react_dir = os.path.join(pwd, "exp_new_reactions")
        #react_dir = os.path.join(pwd, "all_reactions_new")
        # react_dir = os.path.join(pwd, "gascorrected_stable_new_reactions")
        #react_dir = os.path.join(pwd, "revised_reactions_theo_exp")
        react_dir = os.path.join(pwd, "revised_reactions_nobaremetals_2020compat")
        #react_dir = None

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
            react_energy = float(split_line[(parse_i + 2)])
            reaction_split = split_line[(parse_i + 3):]
           
           
            ### converting strings of compounds with compositions ###
            math_symbols = ["+", "->"]
            compounds = [group for group in reaction_split if not (is_number_or_dec(group) or group in math_symbols)]
            #print(f"compounds: {compounds}")
            compounds_compositions = [Composition(cmpd) for cmpd in compounds]
            #print(f"compounds composition: {compounds_compositions}")

            if (Composition(composition_) in compounds_compositions): 
                ### replacing compounds strings in reaction with compositions ###
                reaction_split_temp = copy.deepcopy(reaction_split)
                for i in range(len(compounds)):
                    cmpd = compounds[i] 
                    cmpd_idx = reaction_split.index(cmpd)
                    reaction_split_temp[cmpd_idx] = compounds_compositions[i]

                comp_index = reaction_split_temp.index(Composition(composition_))                
                #print(f"reaction_split_temp: {reaction_split_temp}")

                if (comp_index != 0): comp_coeff = float(reaction_split[(comp_index-1)])
                else: comp_coeff = 1

                comp_all_atoms, comp_nonopen_atoms = get_num_atoms(reaction_split[(comp_index)])
                
                print(f"SE_E_above_hull: {SE_E_above_hull}  comp_all_atoms: {comp_all_atoms}   comp_nonopen_atoms: {comp_nonopen_atoms}  comp_coeff: {comp_coeff}")
                #print(f"correction: {SE_E_above_hull}") #/ (comp_nonopen_atoms * comp_coeff)}")
                print(f"react_energy: {react_energy}")
                react_energy_corrected = react_energy - SE_E_above_hull * comp_all_atoms / comp_nonopen_atoms #/ (comp_nonopen_atoms * comp_coeff)
                print(f"react_energy_corrected: {react_energy_corrected}")
                #react_energy = react_energy_corrected

            if (np.isclose(react_energy, 0, rtol=1e-04)): react_energy = 0


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

 
def get_min_energy_entry(mat_comp_str, entries):
    
    energy_min = sys.maxsize
    entry_min = None    
    mat_comp = Composition(mat_comp_str)
    element_of_interest = mat_comp.elements[0]
    mat_comp_dict = mat_comp.get_el_amt_dict()
    mat_comp_amt = mat_comp_dict[element_of_interest.symbol]
    print(f"element_of_interest: {element_of_interest}")
    print(f"mat_comp_amt: {mat_comp_amt}")

    ### iterating through different entries in MP ###
    for value in entries:    
        if (mat_comp == value.composition):
            value_comp_dict = value.composition.get_el_amt_dict()
            print(f"value_comp_dict: {value_comp_dict}")
            value_amt = value_comp_dict[element_of_interest.symbol]
            print(f"value_amt: {value_amt}")
            scale_factor = mat_comp_amt / value_amt
            print(f"value.composition: {value.composition}")
            print(f"value.energy: {value.energy}")
            energy = value.energy * scale_factor
            print(f"energy postscale: {energy}")

            if (energy < energy_min): 
                energy_min = energy
                entry_min = value

            if (entry_min == None): 
                ValueError(f"MP entry objects being incorrectly accessed or requested material does not exist")
                exit()
    
    return energy_min, entry_min


def compare_baremetals_theoreaction(elements_, relative_mu, composition_, write_dir = None, inputfile=None, comp_benign=None, open_el=None, input_products=None, only_exp_observed=False, _unwanted_elems=[], custom_theo_entries_file = None):
    print(f"compare_baremetals_theoreaction()")
    mpr = MPRester('JxCBb3dqR7jtlimd8RA9DokS7uX8Ru7m')

    ### grabs all entries from materials project that permutations of the ###
    ### provided elements and satisfy the additional criteria ###
    if (only_exp_observed) or (custom_theo_entries_file != None):
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
    no_baremetal_entries = []
    
    if (only_exp_observed):
        for i in range(len(exp_entries)):
            exp_entry = exp_entries[i]
            entry_found = False

            for j in range(len(theo_exp_entries)): 
                theo_exp_entry = theo_exp_entries[j] 
                
                if (exp_entry.reduced_formula == theo_exp_entry.reduced_formula):
                    
                    entries.append(theo_exp_entry)
                    shift_entry = copy.deepcopy(theo_exp_entry)

                    if ((len(_unwanted_elems) != 0) and (len(shift_entry.elements) == 1)):
                        
                        if (shift_entry.elements[0].symbol in unwanted_elements):
                            #print(theo_exp_entry.composition)
                            #print(f"old energy: {theo_exp_entry._energy}")
                            shift_entry._energy += 1000
                            #print(f"new energy: {theo_exp_entry._energy}")

                    no_baremetal_entries.append(shift_entry)
                    entry_found = True
                    
            if not entry_found:
                print(f"ENTRY NOT FOUND: {exp_entry.composition} energy: {exp_entry.energy}")

    else:
        if (custom_theo_entries_file != None): 
            print(f"in reading file")
            for i in range(len(exp_entries)):
                exp_entry = exp_entries[i]
                entry_found = False

                for j in range(len(theo_exp_entries)): 
                    theo_exp_entry = theo_exp_entries[j] 

                    if (exp_entry.reduced_formula == theo_exp_entry.reduced_formula):
                        
                        entries.append(theo_exp_entry)
                        shift_entry = copy.deepcopy(theo_exp_entry)

                        if ((len(_unwanted_elems) != 0) and (len(shift_entry.elements) == 1)):
                            
                            if (shift_entry.elements[0].symbol in unwanted_elements):
                                #print(theo_exp_entry.composition)
                                #print(f"old energy: {theo_exp_entry._energy}")
                                shift_entry._energy += 1000
                                #print(f"new energy: {theo_exp_entry._energy}")

                        no_baremetal_entries.append(shift_entry)
                        entry_found = True
                        
                if not entry_found:
                    print(f"ENTRY NOT FOUND: {exp_entry.composition} energy: {exp_entry.energy}")

            custom_theo_entries = read_select_theoretical_compounds_file(custom_theo_entries_file, elements_)

            for entry in theo_exp_entries:
                if (entry.composition in custom_theo_entries):
                    entries.append(entry)
                    no_baremetal_entries.append(entry)

                else: pass

        else: entries = theo_exp_entries

    baremetal_entries = no_baremetal_entries
    allnontheo_entries = entries
    # theo_entries = theo_exp_entries

    write_dir = "revised_elems_nobaremetals"
    decomp_dir = os.path.join(os.getcwd(), write_dir)
    allatoms_str = '_'.join(elements)
    react_output = allatoms_str 
    stable_output = react_output + "_theo_only_output.txt"
    stable_output = os.path.join(decomp_dir, stable_output)

    if os.path.isdir(decomp_dir): pass
    else: os.mkdir(decomp_dir)
    
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

        CSE_loaded_reactant =  vasprun_obj.get_computed_entry(inc_structure, parameters, data, entry_id)
        #print(f"type(CSE_loaded): {type(CSE_loaded)}")
        baremetal_entries.append(CSE_loaded_reactant)
        allnontheo_entries.append(CSE_loaded_reactant)

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
            
            ### altering energies of gases loaded in ###
            if ("O2_" in product_filename):
                if ("O" in elements_):
                    print(f"adding in O compounds")
                    CSE_loaded._energy = -9.57
                    baremetal_entries.append(CSE_loaded)
                    allnontheo_entries.append(CSE_loaded)

                else: pass

            elif (f"Cl2_" in product_filename):
                if ("Cl" in elements_):
                    print(f"adding in Cl compounds")
                    CSE_loaded._energy = -4.02
                    baremetal_entries.append(CSE_loaded)
                    allnontheo_entries.append(CSE_loaded)
                else: pass
                
            elif (f"H2_" in product_filename):
                if ("H" in elements_):
                    print(f"adding in H compounds")
                    CSE_loaded._energy = -6.99
                    baremetal_entries.append(CSE_loaded)
                    allnontheo_entries.append(CSE_loaded)
                else: pass
                
            ### adding in all other energies ### 
            #else: 
            baremetal_entries.append(CSE_loaded)
            allnontheo_entries.append(CSE_loaded)    
    
    SE_entry = None
    if (inputfile != None):
        SE_entry = CSE_loaded_reactant    
    else:
        for entry in allnontheo_entries:
            if (entry.reduced_formula == Composition(composition_).reduced_formula):
                
                if (SE_entry is None): 
                    SE_entry = entry
                    # print(f"entry.composition: {entry.composition}")
                    baremetal_entries.append(entry)
                    allnontheo_entries.append(entry)

                elif (SE_entry.energy > entry.energy): 
                    SE_entry = entry
                    # print(f"entry.composition: {entry.composition}")
                    baremetal_entries.append(entry)
                    allnontheo_entries.append(entry)


    compat = MaterialsProject2020Compatibility()
    compat.process_entries(baremetal_entries, 
        clean = True,
        verbose = False,
        inplace = True,
        on_error = "ignore")
    compat.process_entries(allnontheo_entries, 
        clean = True,
        verbose = False,
        inplace = True,
        on_error = "ignore")

    
    energy_Li6Zr2O7, entry_Li6Zr2O7 = get_min_energy_entry("Li6Zr2O7", baremetal_entries)
    energy_LLZOtet, entry_LLZOtet = get_min_energy_entry("Li7La3Zr2O12", baremetal_entries)   
    energy_Ta, entry_Ta = get_min_energy_entry("Ta", baremetal_entries)   
    energy_Zr, entry_Zr = get_min_energy_entry("Zr", baremetal_entries)   
    energy_Li2O, entry_Li2O = get_min_energy_entry("Li2O", baremetal_entries)  
    energy_Li5TaO5, entry_Li5TaO5 = get_min_energy_entry("Li5TaO5", baremetal_entries)
    energy_La2Zr2O7, entry_La2Zr2O7 = get_min_energy_entry("La2Zr2O7", baremetal_entries)
    energy_La2O3, entry_La2O3 = get_min_energy_entry("La2O3", baremetal_entries)
    energy_Li, entry_Li = get_min_energy_entry("Li", baremetal_entries)

    filename = inputfile            # xml file name (str)
    ionic_step_skip = None          # read every "ith" frame if set to "i" (int, default = None))
    ionic_step_offset = -1          # start reading at "ith" frame if set to "i" (int, default = 0)
    parse_dos = False               # read in density of states (bool, default = True))
    parse_eigen = False             # read in eigenvalues (bool, default = True))
    parse_projected_eigen = False   # read in projected eigenvals / magnetization (bool, default = True))
    parse_potcar_file = False       # read in potcar (bool, default = True))
    occu_tol = 1e-8                 # minimum tolerance of vbm and cbm (float, default = 1e-8)
    separate_spins= False           # report vbm, cbm, and band gap for each individual spin channel (bool, default = False, must be spin-polarized calc.)
    exception_on_bad_xml = True 
    inc_structure = True          # read in ComputedStructureEntry instead of ComputedEntry (bool, defalut = True)
    parameters = None              # input parameters supported by Vasprun object (list, default = None, default set of params for post-processing given)
    data = None                    # output data to include supported by the Vasprun object (dict, default = None)
    entry_id = None      

    filename_LLTZO = "vasprun_LLTZO.xml"   # xml file name (str)
    vasprun_obj_LLZTO = Vasprun(filename_LLTZO, ionic_step_skip, ionic_step_offset, parse_dos, parse_eigen, parse_projected_eigen, 
                    parse_potcar_file, occu_tol, separate_spins, exception_on_bad_xml)
    entry_LLTZO =  vasprun_obj_LLZTO.get_computed_entry(inc_structure, parameters, data, entry_id)
    print(f"entry_LLTZO.energy: {entry_LLTZO.energy}")
    LLTZO_dict = entry_LLTZO.composition.get_el_amt_dict()
    value_amt = LLTZO_dict["Ta"]
    print(f"value_amt: {value_amt}")
    scale_factor = 1 / value_amt
    energy_LLTZO = entry_LLTZO.energy * scale_factor

    filename_LiLa2TaO6 = "vasprun_LiLa2TaO6.xml"   # xml file name (str)
    vasprun_obj_LiLa2TaO6 = Vasprun(filename_LiLa2TaO6, ionic_step_skip, ionic_step_offset, parse_dos, parse_eigen, parse_projected_eigen, 
                    parse_potcar_file, occu_tol, separate_spins, exception_on_bad_xml)
    entry_LiLa2TaO6 = vasprun_obj_LiLa2TaO6.get_computed_entry(inc_structure, parameters, data, entry_id)
    print(f"entry_LiLa2TaO6.energy: {entry_LiLa2TaO6.energy}")
    LiLa2TaO6_dict = entry_LiLa2TaO6.composition.get_el_amt_dict()
    value_amt = LiLa2TaO6_dict["Li"]
    print(f"value_amt: {value_amt}")
    scale_factor = 1 / value_amt
    energy_LiLa2TaO6 = entry_LiLa2TaO6.energy * scale_factor

    baremetal_pd = PhaseDiagram(baremetal_entries)
    theo_exp_pd = PhaseDiagram(allnontheo_entries)

    print(f"\n")
    
    # Get the chemical potential of the pure subtance.
    mu = theo_exp_pd.get_transition_chempots(Element(open_el))[0]
    
    # Set the chemical potential in the elemental reservoir.
    chempots = {open_el: relative_mu + mu}

    # Build the grand potential phase diagram
    theo_exp_gpd = GrandPotentialPhaseDiagram(entries=allnontheo_entries, chempots=chempots) #, elements=elements_)
    baremetal_gpd = GrandPotentialPhaseDiagram(entries=baremetal_entries, chempots=chempots) #, elements=elements_)
    
    # plotter_theoexp = PDPlotter(theo_exp_pd).get_plot()
    # plotter_theoexp.show()
    # plotter_baremetals = PDPlotter(baremetal_pd).get_plot()
    # plotter_baremetals.show()

    comp_SE = Composition(composition_)
    theo_exp_reacs = GrandPotentialInterfacialReactivity(c1=comp_SE, c2=comp_benign, grand_pd=theo_exp_gpd, pd_non_grand=theo_exp_pd, norm=True, include_no_mixing_energy=False, use_hull_energy=False)
    baremetal_reacs = GrandPotentialInterfacialReactivity(c1=comp_SE, c2=comp_benign, grand_pd=baremetal_gpd, pd_non_grand=baremetal_pd, norm=True, include_no_mixing_energy=False, use_hull_energy=False)
    
    
    print(f"theoretical reaction")
    critical_rxns = [
        OrderedDict([
            ("Atomic fraction", round(ratio, 4)),
            ("Reaction equation", rxn),
            ("E$_{rxt}$ per mol equation (kJ/mol)", round(rxn_energy, 4)),
            ("E$_{rxt}$ per reactant atom (eV/atom)", round(reactivity, 4)),
        ])

        for _, ratio, reactivity, rxn, rxn_energy in theo_exp_reacs.get_kinks()]
    GPIR_interface_reaction_table = pandas.DataFrame(critical_rxns)
    with pandas.option_context('display.max_rows', None, 'display.max_columns', None, 'display.max_colwidth', None):  # more options can be specified also
        print(f"GPIR_interface_reaction_table: {GPIR_interface_reaction_table}")
    
    print(f"bare metal reaction")
    critical_rxns = [
        OrderedDict([
            ("Atomic fraction", round(ratio, 4)),
            ("Reaction equation", rxn),
            ("E$_{rxt}$ per mol equation (kJ/mol)", round(rxn_energy, 4)),
            ("E$_{rxt}$ per reactant atom (eV/atom)", round(reactivity, 4)),
        ])

        for _, ratio, reactivity, rxn, rxn_energy in baremetal_reacs.get_kinks()]
    GPIR_interface_reaction_table = pandas.DataFrame(critical_rxns)
    with pandas.option_context('display.max_rows', None, 'display.max_columns', None, 'display.max_colwidth', None):  # more options can be specified also
        print(f"GPIR_interface_reaction_table: {GPIR_interface_reaction_table}")
    

    '''
    ### calculate reaction energies ###
    
    ### cole equation ###
    deltaE_cole = ( 11*energy_LLZOtet + 5*energy_Li2O + energy_Li6Zr2O7 + 6*energy_LiLa2TaO6) - ( 15 * energy_LLTZO )

    ### christian equation ###
    deltaE_christian = ( 3*energy_La2O3 + 3*energy_Zr + 15*energy_Li2O + energy_Ta) - ( 17*energy_Li + energy_LLTZO )

    ### miara equation ###
    deltaE_miara = ( 2*energy_La2Zr2O7 + 34*energy_La2O3 + 16*energy_Li6Zr2O7 + 12*energy_Li5TaO5) - ( 24 * energy_LLTZO )

    print(f"deltaE_cole: {deltaE_cole}")
    print(f"deltaE_christian: {deltaE_christian / 44}")
    print(f"deltaE_miara: {deltaE_miara}")
    print(f"deltaE_miara / num_atoms: {deltaE_miara / 47}") 
    '''

    if (write_dir != None):
        allatoms_str = '_'.join(elements)
        react_output = allatoms_str + f"_{round(relative_mu, 3)}volts_reactions.txt"

        pwd = os.getcwd()
        react_dir = os.path.join(pwd, "compare_theo_baremetal_reacenergy")

        if os.path.isdir(react_dir): pass
        else: os.mkdir(react_dir)
        react_output = os.path.join(react_dir, react_output)

        print(f"theoretical reaction")
        react_file = open(react_output, "w")
        react_file.write(f"{composition_} theoretical decomposition reactions\n--------\n")

        theo_reac_energies = []
        nobaremetal_reac_energies = []

        for _, value in theo_exp_reacs.labels.items():
            split_line = value.split()
            parse_i = 0
            while ("eV/atom" not in split_line[parse_i]): parse_i += 1
            reac_energy_printed = float(split_line[(parse_i + 2)])
            reaction_chunk = " ".join(split_line[(parse_i + 3):])
            
            reactant_dict, product_dict, _ = parse_reaction(reaction_chunk)

            reactant_energy = 0
            reactant_grand_pot = 0
            reac_atoms_sum = 0
            reac_atoms_nonopen_sum = 0

            product_energy = 0
            product_grand_pot = 0
            product_atoms_sum = 0
            product_atoms_nonopen_sum = 0

            ### finding energy contributions of reactants
            for compound, coeff in reactant_dict.items():
                compound_noopen = Composition({k: v for k, v in Composition(compound).items() if k not in [Element("Na"), Element("Li")]})
                # entry_found = identify_lowest_e_entry(compound_noopen, theo_exp_entries)
                # energy = float(entry_found[0]) * coeff

                # normalized_energy = self._get_energy(mixing_ratio)
                #if (compound not in ["Li1", "Na1"]): energy = theo_exp_gpd.get_hull_energy(compound_noopen) * coeff
                #else: energy = 0
                # Gets balanced reaction at kinks
                # rxt_energy = normalized_energy * self._get_elem_amt_in_rxn(rxt) / n_atoms
                num_atoms, num_non_open_atoms = get_num_atoms(compound)

                if (compound not in ["Li1", "Na1"]): energy = theo_exp_reacs._get_grand_potential(Composition(compound)) 
                else: energy = 0

                reactant_energy += energy #/ num_atoms
                #if (compound not in ["Li1", "Na1"]): print(f"theo_exp_gpd.get_hull_energy(compound_noopen): {energy / (coeff * num_non_open_atoms)}")
                #else: print(f"theo_exp_gpd.get_hull_energy: 0")
                #print(f"energy: {energy}")
                #print(f"compound: {compound}   compound_noopen: {compound_noopen}  num_atoms: {num_atoms}   num_non_open_atoms: {num_non_open_atoms}")
                
                if (compound not in ["Li1", "Na1"]): reactant_grand_pot += energy # * num_non_open_atoms / num_atoms
                reac_atoms_sum += num_atoms * coeff 
                reac_atoms_nonopen_sum += num_non_open_atoms * coeff 
            

            ### finding energy contributions of products
            for compound, coeff in product_dict.items():
                
                num_atoms, num_non_open_atoms = get_num_atoms(compound)
                compound_noopen = Composition({k: v for k, v in Composition(compound).items() if k not in [Element("Na"), Element("Li")]})

                if (compound not in ["Li1", "Na1"]): energy = theo_exp_reacs._get_grand_potential(Composition(compound)) * coeff 
                else: energy = 0
                
                # print(f"energy: {energy}")
                
                product_energy += energy 
                
                if (compound not in ["Li1", "Na1"]): product_grand_pot += energy  * int(num_non_open_atoms) / int(reac_atoms_nonopen_sum)#
                product_atoms_sum += num_atoms * coeff
                product_atoms_nonopen_sum += num_non_open_atoms * coeff
                

            #print(f"product_grand_pot: {product_grand_pot}   reactant_grand_pot: {reactant_grand_pot}")
            reac_energy = (product_energy - reactant_energy) 
            reac_grand_pot = (product_grand_pot - reactant_grand_pot) 

            split_value = split_line[(parse_i + 3):]
            join_value = " ".join(split_value) + f" normal: {reac_energy_printed} eV/atom,  calc E: {reac_energy},  calc GP: {reac_grand_pot}"
            print(f"write: {join_value}")
            react_file.write(f"{join_value}\n")
            theo_reac_energies.append(tuple([reac_energy, reac_grand_pot]))

        print(f"bare metal reaction")
        react_file.write(f"\n")
        react_file.write(f"{composition_} no bare metal decomposition reactions\n--------\n")
        
        for _, value in baremetal_reacs.labels.items():
            split_line = value.split()
            parse_i = 0
            while ("eV/atom" not in split_line[parse_i]): parse_i += 1
            reac_energy_printed = float(split_line[(parse_i + 2)])
            reaction_chunk = " ".join(split_line[(parse_i + 3):])
            
            reactant_dict, product_dict, _ = parse_reaction(reaction_chunk)

            reactant_energy = 0
            reactant_grand_pot = 0
            reac_atoms_sum = 0
            reac_atoms_nonopen_sum = 0

            product_energy = 0
            product_grand_pot = 0
            product_atoms_sum = 0
            product_atoms_nonopen_sum = 0

            ### finding energy contributions of reactants
            for compound, coeff in reactant_dict.items():
                compound_noopen = Composition({k: v for k, v in Composition(compound).items() if k not in [Element("Na"), Element("Li")]})
                
                # Gets balanced reaction at kinks
                # rxt_energy = normalized_energy * self._get_elem_amt_in_rxn(rxt) / n_atoms
                num_atoms, num_non_open_atoms = get_num_atoms(compound)

                if (compound not in ["Li1", "Na1"]): energy = theo_exp_reacs._get_grand_potential(Composition(compound)) 
                else: energy = 0

                reactant_energy += energy #/ num_atoms
                
                if (compound not in ["Li1", "Na1"]): reactant_grand_pot += energy # * num_non_open_atoms / num_atoms
                reac_atoms_sum += num_atoms * coeff 
                reac_atoms_nonopen_sum += num_non_open_atoms * coeff 
            

            ### finding energy contributions of products
            for compound, coeff in product_dict.items():
                
                num_atoms, num_non_open_atoms = get_num_atoms(compound)
                compound_noopen = Composition({k: v for k, v in Composition(compound).items() if k not in [Element("Na"), Element("Li")]})

                if (compound not in ["Li1", "Na1"]): energy = theo_exp_reacs._get_grand_potential(Composition(compound)) * coeff 
                else: energy = 0
                 
                product_energy += energy 
                
                if (compound not in ["Li1", "Na1"]): product_grand_pot += energy  * int(num_non_open_atoms) / int(reac_atoms_nonopen_sum)#
                product_atoms_sum += num_atoms * coeff
                product_atoms_nonopen_sum += num_non_open_atoms * coeff
                

            #print(f"product_grand_pot: {product_grand_pot}   reactant_grand_pot: {reactant_grand_pot}")
            reac_energy = (product_energy - reactant_energy) 
            reac_grand_pot = (product_grand_pot - reactant_grand_pot) 

            split_value = split_line[(parse_i + 3):]
            join_value = " ".join(split_value) + f" normal: {reac_energy_printed} eV/atom,  calc E: {reac_energy},  calc GP: {reac_grand_pot}"
            print(f"write: {join_value}")
            react_file.write(f"{join_value}\n")
            nobaremetal_reac_energies.append(tuple([reac_energy, reac_grand_pot]))

        react_file.close()
        print("\n\n\n")

    return theo_reac_energies[-1][1], nobaremetal_reac_energies[-1][1]


mu_vals = np.arange(-1, 1.01, 0.1) 
#mu_vals = [0]
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

plot_dirpath = "comparison_plots"
energy_difference = np.zeros((len(mu_vals),len(electrolyte_compounds),2))
idx = 0

for mu_in in mu_vals:

    print(f"mu_in {mu_in} volts")
    input_filename = "vasprun_new_LLZTO.xml" 
    inputproducts = ["vasprun_new_LiLa2TaO6.xml", "vasprun_O2.xml", "vasprun_Cl2.xml", "vasprun_H2.xml"]
    '''
    elements = ['Li', 'La', 'Zr', 'O', 'Ta'] # LLTZO
    compare_baremetals_theoreaction(elements, mu_in, "Li13La6Zr3TaO24", comp_benign=Composition("Li2O"), open_el="Li", 
         write_dir = outputdir, only_exp_observed=only_exp_observed_val, inputfile = input_filename, input_products=inputproducts, _unwanted_elems=unwanted_elements, custom_theo_entries_file=custom_theoretical_compounds)
    
    elements = ['Li', 'La', 'Zr', 'O', 'Ta'] # LLTZO
    theo_E, baremetal_E = compare_baremetals_theoreaction(elements, mu_in, "Li13La6Zr3TaO24", comp_benign=Composition("Li2O"), open_el="Li", 
         write_dir = outputdir, only_exp_observed=only_exp_observed_val, inputfile = input_filename, input_products=inputproducts, _unwanted_elems=unwanted_elements, custom_theo_entries_file=custom_theoretical_compounds)
    '''
    

    # inputproducts = ["vasprun_O2.xml", "vasprun_Cl2.xml", "H2_vasprun.xml"]
    input_filename = None
    inputproducts = ["vasprun_O2.xml"]
    elements = ['Li', 'La', 'Zr', 'O'] # LLZO    
    theo_E, baremetal_E = compare_baremetals_theoreaction(elements, mu_in, "Li7La3Zr2O12", comp_benign=Composition("Li2O"), open_el="Li", 
        write_dir = outputdir, only_exp_observed=only_exp_observed_val, inputfile = input_filename, input_products=inputproducts, _unwanted_elems=unwanted_elements, custom_theo_entries_file=custom_theoretical_compounds)

    inputproducts = None
    elements = ['Li', 'Ge', 'P', 'S'] # Li10GeP2S12
    theo_E, baremetal_E = compare_baremetals_theoreaction(elements, mu_in, "Li10GeP2S12", comp_benign=Composition("Li2S"), open_el="Li", 
        write_dir = outputdir, only_exp_observed=only_exp_observed_val, inputfile = input_filename, input_products=inputproducts, _unwanted_elems=unwanted_elements, custom_theo_entries_file=custom_theoretical_compounds)

    
    elements = ['Na', 'Sb', 'S'] # Na3SbS4
    theo_E, baremetal_E = compare_baremetals_theoreaction(elements, mu_in, "Na3SbS4", comp_benign=Composition("Na2S"), open_el="Na", 
        write_dir = outputdir, only_exp_observed=only_exp_observed_val, inputfile = input_filename, input_products=inputproducts, _unwanted_elems=unwanted_elements, custom_theo_entries_file=custom_theoretical_compounds)

    
    inputproducts = ["vasprun_Cl2.xml"]
    elements = ['Li', 'In', 'Cl'] # Li3InCl6
    theo_E, baremetal_E = compare_baremetals_theoreaction(elements, mu_in, 'Li3InCl6', comp_benign=Composition("LiCl"), open_el="Li", 
        write_dir = outputdir, only_exp_observed=only_exp_observed_val, inputfile = input_filename, input_products=inputproducts, _unwanted_elems=unwanted_elements, custom_theo_entries_file=custom_theoretical_compounds)


    inputproducts = ["vasprun_O2.xml"]
    elements = ['Na', 'O', 'Br'] # Na3OBr 
    theo_E, baremetal_E = compare_baremetals_theoreaction(elements, mu_in, "Na3OBr", comp_benign=Composition("Na2O"), open_el="Na", 
        write_dir = outputdir, only_exp_observed=only_exp_observed_val, inputfile = input_filename, input_products=inputproducts, _unwanted_elems=unwanted_elements, custom_theo_entries_file=custom_theoretical_compounds)

    
    inputproducts = ["vasprun_Cl2.xml"]
    elements = ['Li', 'P', 'S', 'Cl'] # Li6P1S5Cl1 
    theo_E, baremetal_E = compare_baremetals_theoreaction(elements, mu_in, "Li6P1S5Cl1", comp_benign=Composition("Li2S"), open_el="Li", 
        write_dir = outputdir, only_exp_observed=only_exp_observed_val, inputfile = input_filename, input_products=inputproducts, _unwanted_elems=unwanted_elements, custom_theo_entries_file=custom_theoretical_compounds)


    inputproducts = None
    elements = ['Li', 'P', 'S'] # Li7P3S11 
    theo_E, baremetal_E = compare_baremetals_theoreaction(elements, mu_in, "Li7P3S11", comp_benign=Composition("Li2S"), open_el="Li", 
        write_dir = outputdir, only_exp_observed=only_exp_observed_val, inputfile = input_filename, input_products=inputproducts, _unwanted_elems=unwanted_elements, custom_theo_entries_file=custom_theoretical_compounds)


    inputproducts = ["vasprun_O2.xml"]
    elements = ['Na', 'Al', 'O'] ### na-beta composition NaAl11O17 
    theo_E, baremetal_E = compare_baremetals_theoreaction(elements, mu_in, "NaAl11O17", comp_benign=Composition("Na2O"), open_el="Na", 
        write_dir = outputdir, only_exp_observed=only_exp_observed_val, inputfile = input_filename, input_products=inputproducts, _unwanted_elems=unwanted_elements, custom_theo_entries_file=custom_theoretical_compounds)


    inputproducts = ["vasprun_Cl2.xml", "vasprun_O2.xml", "vasprun_H2.xml"]
    elements = ['Li', 'O', 'H', 'Cl'] ### Lithium antiperovskite composition Li2OHCl 
    theo_E, baremetal_E = compare_baremetals_theoreaction(elements, mu_in, "Li2OHCl", comp_benign=Composition("Li2O"), open_el="Li", 
        write_dir = outputdir, only_exp_observed=only_exp_observed_val, inputfile = input_filename, input_products=inputproducts, _unwanted_elems=unwanted_elements, custom_theo_entries_file=custom_theoretical_compounds)

    inputproducts = ["vasprun_O2.xml"]
    elements = ['Na', 'Zr', 'Si', 'O'] ### Lithium antiperovskite composition Li2OHCl 
    theo_E, baremetal_E = compare_baremetals_theoreaction(elements, mu_in, "Na4Zr2Si3O12", comp_benign=Composition("Na2O"), open_el="Na", 
        write_dir = outputdir, only_exp_observed=only_exp_observed_val, inputfile = input_filename, input_products=inputproducts, _unwanted_elems=unwanted_elements, custom_theo_entries_file=custom_theoretical_compounds)
