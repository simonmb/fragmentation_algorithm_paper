# -*- coding: utf-8 -*-
''' Code to execute and evaluate fragmenting molecules into molecular subgroups

MIT License

Copyright (C) 2019, Simon Mueller <simon.mueller@tuhh.de>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.'''

import operator
from fragmenter import fragmenter
try:
    from rdkit import Chem
except:
    raise Exception('RDkit has to be installed.')

def fix_SMILES(SMILES):
    SMILES = SMILES.replace('[1H]', '')
    SMILES = SMILES.replace('[2H]', '')
    SMILES = SMILES.replace('[3H]', '')
    SMILES = SMILES.replace('[H]', '')
    return SMILES

def info_to_CSV(inchikey, SMILES, pubchem_id, fragmentation):
    
    fragmentation_array = []
    for group_number, amount in fragmentation.items():
        fragmentation_array.append(str(group_number) + ":" + str(amount))
    
    return inchikey + "," + SMILES + "," + pubchem_id + "," + "|".join(fragmentation_array)

def CSV_to_info(CSV_line, has_fragmentation = False):
    CSV_line = CSV_line.replace('\n', '')
    array = CSV_line.split(',')
    
    fragmentation = {}
    
    if has_fragmentation:
        fragmentation_array = array[3].split('|')
        for match_str in fragmentation_array:
            array2 = match_str.split(':')
            group_number = int(array2[0])
            amount = int(array2[1])
            
            fragmentation[group_number] = amount
        
    return array[0], array[1], array[2], fragmentation
    

def function_to_choose_fragmentation(fragmentations):
    fragmentations_descriptors = {}
    i = 0
    for fragmentation in fragmentations:
        fragmentations_descriptors[i] = [len(fragmentation)]
        i += 1
    
    sorted_fragmentations_dict = sorted(fragmentations_descriptors.items(), key=operator.itemgetter(1))

    return fragmentations[sorted_fragmentations_dict[0][0]]

def is_fragmentation_equal_to_other_fragmentation(fragmentation, other_fragmentation):
    
    for group_number, amount in fragmentation.items():
        if other_fragmentation.has_key(group_number):
            if fragmentation[group_number] != other_fragmentation[group_number]:
                return False
    return True

def log_structure_results(f, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB, status = ''):
    
    f.write('https://pubchem.ncbi.nlm.nih.gov/compound/' + pubchem_id + '#section=2D-Structure\n')
    f.write(SMILES + '\n')
    f.write(inchikey + '\n')
    f.write('\n' + 'Fragmentation was successfull: ' + str(success) + '\n')
    
    if status != '':
        f.write(status + '\n')
    
    if success:
        f.write('Fragmentation from the algorithm:\n')
        sorted_group_number = sorted(fragmentation.iterkeys())
        
        for group_number in sorted_group_number:
            f.write((UNIFAC_SMARTS[group_number - 1][0]).ljust(12, ' ') + '\t' + str(group_number).ljust(8, ' ') + str(fragmentation[group_number]).ljust(8, ' ') + '\n')
    
    f.write('\n')
    
    if len(fragmentation_reference_DB) > 0:
        f.write('Fragmentation from the reference database:\n')
        sorted_group_number = sorted(fragmentation_reference_DB.iterkeys())
        
        for group_number in sorted_group_number:
            f.write((UNIFAC_SMARTS[group_number - 1][0]).ljust(12, ' ') + '\t' + str(group_number).ljust(8, ' ') + str(fragmentation_reference_DB[group_number]).ljust(8, ' ') + '\n')
    
        
    f.write('\n\n')   


UNIFAC_SMARTS =  [
("CH3", "[CH3;X4]"),
("CH2", "[CH2;X4]"),
("CH", "[CH1;X4]"),
("C", ["[CH0;X4]", "[CH0;X3]"]),
("CH2=CH", "[CH2]=[CH]"),
("CH=CH", "[CH]=[CH]"),
("CH2=C", ["[CH2]=[C]", "[CH2]=[c]"]),
("CH=C", ["[CH]=[CH0]", "[CH]=[cH0]"]),
("ACH", "[cH]"),
("AC", "[cH0]"),
("ACCH3", "[c][CH3;X4]"),
("ACCH2", "[c][CH2;X4]"),
("ACCH", "[c][CH;X4]"),
('OH', "[OH]"),
('CH3OH', "[CH3][OH]"),
('H2O', "[OH2]"),
('ACOH', "[c][OH]"),
("CH3CO", "[CH3][CH0]=O"),
("CH2CO", "[CH2][CH0]=O"),
("CH=O", "[CH]=O"),
("CH3COO", "[CH3]C(=O)[OH0]"),
("CH2COO", "[CH2]C(=O)[OH0]"),
("HCOO", "[CH](=O)[OH0]"),
("CH3O", "[CH3][OH0]"),
("CH2O", "[CH2][OH0]"),
("CHO", "[CH][OH0]"),
("THF", "[CH2;R][OH0]"),
("CH3NH2", "[CH3][NH2]"),
("CH2NH2", "[CH2][NH2]"),
("CHNH2", "[CH][NH2]"),
("CH3NH", "[CH3][NH]"),
("CH2NH", "[CH2][NH]"),
("CHNH", "[CH][NH]"),
("CH3N",  ["[CH3][N]", "[CH3][n]"]),
("CH2N", "[CH2][N]"),
("ACNH2", "[c][NH2]"),
("C5H5N", "n1[cH][cH][cH][cH][cH]1"),
("C5H4N", ["n1[c][cH][cH][cH][cH]1",
           "n1[cH][c][cH][cH][cH]1",
           "n1[cH][cH][c][cH][cH]1"]),
('C5H3N', ["n1[c][c][cH][cH][cH]1",
           "n1[c][cH][c][cH][cH]1",
           "n1[c][cH][cH][c][cH]1",
           "n1[c][cH][cH][cH][c]1",
           "n1[cH][c][c][cH][cH]1",
           "n1[cH][c][cH][c][cH]1"]),
("CH3CN", "[CH3]C#N"),
("CH2CN", "[CH2]C#N"),
("COOH", "C(=O)[OH]"),
("HCOOH", "[CH](=O)[OH]"),
("CH2Cl", "[CH2]Cl"),
("CHCl", "[CH]Cl"),
("CCl", "[CH0]Cl"),
("CH2Cl2", "[CH2](Cl)Cl"),
("CHCl2", "[CH](Cl)Cl"),
("CCl2", "C(Cl)Cl"),
("CHCl3", "[CH](Cl)(Cl)Cl"),
("CCl3", "C(Cl)(Cl)(Cl)"),
("CCl4", "C(Cl)(Cl)(Cl)(Cl)"),
("ACCl", "[c]Cl"),
("CH3NO2", "[CH3][N+](=O)[O-]"),
("CH2NO2", "[CH2][N+](=O)[O-]"),
("CHNO2", "[CH][N+](=O)[O-]"),
("ACNO2", "[c][N+](=O)[O-]"),
("CS2", "C(=S)=S"),
("CH3SH", "[CH3][SH]"),
("CH2SH", "[CH2][SH]"),
("Furfural", "O=[CH]c1[cH][cH][cH]o1"),
("DOH", "[OH][CH2][CH2][OH]"),
("I", "[IH0]"),
("Br", "[BrH0]"),
("CH#C", "[CH]#C"),
("C#C", "C#C"),
("DMSO", "[CH3]S(=O)[CH3]"),
("ACRY", "[CH2]=[CH1][C]#N"),
("Cl(C=C)", "[$(Cl[C]=[C])]"),
("C=C", "[CH0]=[CH0]"),
("ACF", "[c]F"),
("DMF", "[CH](=O)N([CH3])[CH3]"),
("HCON(CH2)2", ["[CH](=O)N([CH2])[CH2]", "[CH](=O)N([CH2])[CH3]"]),
("CF3", "C(F)(F)F"),
("CF2", "C(F)F"),
("CF", "[C]F"),
("COO", ["[CH0](=O)[OH0]", "[cH0](=O)[oH0]"]),  
("SiH3", "[SiH3]"),
("SiH2", "[SiH2]"),
("SiH", "[SiH]"),
("Si", "[Si]"),
("SiH2O", "[SiH2][OH0]"),
("SiHO", "[SiH][OH0]"),
("SiO", "[Si][OH0]"),
("NMP", "[CH3]N1[CH2][CH2][CH2]C(=O)1"),
("CCl3F", "C(Cl)(Cl)(Cl)F"),
("CCl2F", "C(Cl)(Cl)F"),
("HCCl2F", "[CH](Cl)(Cl)F"),
("HCClF", "[CH](Cl)F"),
("CClF2", "C(Cl)(F)F"),
("HCClF2", "[CH](Cl)(F)F"),
("CClF3", "C(Cl)(F)(F)F"),
("CCl2F2", "C(Cl)(Cl)(F)F"),
("CONH2", "C(=O)[NH2]"),
("CONHCH3", "C(=O)[NH][CH3]"),
("CONHCH2", "C(=O)[NH][CH2]"),
("CON(CH3)2", "C(=O)N([CH3])[CH3]"),
("CONCH3CH2", "C(=O)N([CH3])[CH2]"),
("CON(CH2)2", "C(=O)N([CH2])[CH2]"),
("C2H5O2", "[OH0;!$(OC=O);!R][CH2;!R][CH2;!R][OH]"),
("C2H4O2", ["[OH0;!$(OC=O);!R][CH;!R][CH2;!R][OH]",
              "[OH0;!$(OC=O);!R][CH2;!R][CH;!R][OH]"]),
("CH3S", "[CH3]S"),
("CH2S", "[CH2]S"),
("CHS", "[CH]S"),
("MORPH", "[CH2]1[CH2][NH][CH2][CH2]O1"),
("C4H4S", "[cH]1[cH][s;X2][cH][cH]1"),
('C4H3S', ["[c]1[cH][s;X2][cH][cH]1",
           "[cH]1[c][s;X2][cH][cH]1"]),
('C4H2S', ["[c]1[c][s;X2][cH][cH]1",
           "[c]1[cH][s;X2][cH][c]1",
           "[cH]1[c][s;X2][c][cH]1",
           "[cH]1[c][s;X2][cH][c]1"]),
("NCO", "N=C=O"),
("H2COCH", ""),    # "[CH2]1[CH]O1"
("HCOCH", ""),     # "[CH]1[CH]O1"
("COCH", ""),      # "C1[CH]O1"
("H2COCH2", ""),   # "[CH2]1[CH2]O1"
("OCOCO", ""),     # "C(=O)OC(=O)"),
("(CH3O)2CO", ""), # "[CH3]OC(=O)O[CH3]"
("(CH2O)2CO", ""), # "[CH2]OC(=O)O[CH2]"
("CH3OCH2OCO", ""),# "[CH3]OC(=O)O[CH2]"
("(CH2)2SU", "[CH2]S(=O)(=O)[CH2]"),
("CH2CHSU", "[CH2]S(=O)(=O)[CH]")]

# get the fragmentation scheme in the format necessary
fragmentation_scheme = {i+1: j[1] for i, j in enumerate(UNIFAC_SMARTS)}

# sort the fragmentation scheme according to the descriptors
pattern_descriptors = {}
for group_number, SMARTS in fragmentation_scheme.items():
    if type(SMARTS) is list:
        SMARTS = SMARTS[0]
    
    if SMARTS != "":
        pattern = fragmenter.get_mol_with_properties_from_SMARTS(SMARTS)
        
        pattern_descriptors[group_number] = [pattern.GetUnsignedProp('n_available_bonds') == 0, \
                                  (pattern.GetBoolProp('is_simple_atom_on_c') or pattern.GetBoolProp('is_simple_atom')), \
                                  pattern.GetUnsignedProp('n_atoms_defining_SMARTS'),  
                                  pattern.GetUnsignedProp('n_available_bonds') == 1, \
                                  fragmenter.get_heavy_atom_count(pattern) - pattern.GetUnsignedProp('n_carbons'), \
                                  pattern.GetBoolProp('has_atoms_in_ring'), \
                                  pattern.GetUnsignedProp('n_triple_bonds'), \
                                  pattern.GetUnsignedProp('n_double_bonds')]

sorted_pattern_descriptors = sorted(pattern_descriptors.items(), key=operator.itemgetter(1), reverse=True)
sorted_group_numbers = [i[0] for i in sorted_pattern_descriptors]


# first step: fragment reference database and compare with the results
reference_DB = []
with open('reference_DB.csv') as f:
    for line in f.readlines():
        reference_DB.append(CSV_to_info(line, True))

reference_DB_fragmentation_stats = {}

simple_fragmenter_fragmented = []
simple_fragmenter_fragmented_and_equal_to_reference_DB = []
complete_fragmenter_fragmented = []
complete_fragmenter_fragmented_and_equal_to_reference_DB = []

right_size_for_complete_fragmenter = []

simple_fragmenter = fragmenter(fragmentation_scheme, 'simple')
complete_fragmenter = fragmenter(fragmentation_scheme, 'complete', 20, function_to_choose_fragmentation)

# without sorting the patterns
print '####################################################################'
print 'Fragmenting the reference database without the patterns sorted (simple and complete algorithm)'
i_structure = 0
f_simple = open('reference_DB_simple_fragmentation_without_pattern_sorting_results.log','w+')
f_complete = open('reference_DB_complete_fragmentation_without_pattern_sorting_results.log','w+')
for inchikey, SMILES, pubchem_id, fragmentation_reference_DB in reference_DB:
    
    i_structure = i_structure + 1
    if i_structure % 2000 == 0:
        print '{:2.1f} .'.format((100.0 * i_structure) / len(reference_DB)),
        
    lines = []
    
    
    for group_number, amount in fragmentation_reference_DB.items():
        if not reference_DB_fragmentation_stats.has_key(group_number):
            reference_DB_fragmentation_stats[group_number] = 0
            
        reference_DB_fragmentation_stats[group_number] += amount
    
    fragmentation, success = simple_fragmenter.fragment(SMILES)
    if success:
        simple_fragmenter_fragmented.append(inchikey)
        if is_fragmentation_equal_to_other_fragmentation(fragmentation, fragmentation_reference_DB):
            simple_fragmenter_fragmented_and_equal_to_reference_DB.append(inchikey)
    
    log_structure_results(f_simple, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB)
    
    n_heavy_atoms  = 0
    for sub_SMILES in SMILES.split("."):
        n_heavy_atoms = max(n_heavy_atoms, fragmenter.get_heavy_atom_count(Chem.MolFromSmiles(sub_SMILES)))
    
    if n_heavy_atoms <= 20:
        right_size_for_complete_fragmenter.append(inchikey)
        fragmentation, success = complete_fragmenter.fragment(SMILES)
        if success:
            complete_fragmenter_fragmented.append(inchikey)
            if is_fragmentation_equal_to_other_fragmentation(fragmentation, fragmentation_reference_DB):
                complete_fragmenter_fragmented_and_equal_to_reference_DB.append(inchikey)
        
        log_structure_results(f_complete, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB)
    else:
        log_structure_results(f_complete, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB, 'Structure was skipped because it is larger than 20 atoms.')


f_simple.close()
f_complete.close()


print ''
print 'N_structures(simple): ' + str(len(reference_DB))
print 'N_fragmented(simple): ' + str(len(simple_fragmenter_fragmented)) + "(" + str((1.0 * len(simple_fragmenter_fragmented)) / len(reference_DB)) + ")"
print 'N_fragmented_and_equal(simple): ' + str(len(simple_fragmenter_fragmented_and_equal_to_reference_DB)) + "(" + str((1.0 * len(simple_fragmenter_fragmented_and_equal_to_reference_DB)) / len(reference_DB)) + ")"
print ''
print 'N_structures(complete):' + str(len(right_size_for_complete_fragmenter))
print 'N_fragmented(complete): ' + str(len(complete_fragmenter_fragmented)) + "(" + str((1.0 * len(complete_fragmenter_fragmented)) / len(right_size_for_complete_fragmenter)) + ")"
print 'N_fragmented_and_equal(complete): ' + str(len(complete_fragmenter_fragmented_and_equal_to_reference_DB)) + "(" + str((1.0 * len(complete_fragmenter_fragmented_and_equal_to_reference_DB)) / len(right_size_for_complete_fragmenter)) + ")"
print ''
print ''

# with sorting the patterns
simple_fragmenter.fragmentation_scheme_order = sorted_group_numbers
complete_fragmenter.fragmentation_scheme_order = sorted_group_numbers

simple_fragmenter_sorted_fragmented = []
simple_fragmenter_sorted_fragmented_and_equal_to_reference_DB = []
complete_fragmenter_sorted_fragmented = []
complete_fragmenter_sorted_fragmented_and_equal_to_reference_DB = []

right_size_for_complete_fragmenter2 = []
print '####################################################################'
print 'Fragmenting the reference database with the patterns sorted (simple and complete algorithm)'
i_structure = 0
f_simple = open('reference_DB_simple_fragmentation_with_pattern_sorting_results.log','w+')
f_complete = open('reference_DB_complete_fragmentation_with_pattern_sorting_results.log','w+')
for inchikey, SMILES, pubchem_id, fragmentation_reference_DB in reference_DB:   
             
    i_structure = i_structure + 1
    if i_structure % 2000 == 0:
        print '{:2.1f} .'.format((100.0 * i_structure) / len(reference_DB)),
        
    fragmentation, success = simple_fragmenter.fragment(SMILES)
    if success:
        simple_fragmenter_sorted_fragmented.append(inchikey)
        if is_fragmentation_equal_to_other_fragmentation(fragmentation, fragmentation_reference_DB):
            simple_fragmenter_sorted_fragmented_and_equal_to_reference_DB.append(inchikey)
            
    log_structure_results(f_simple, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB)
    
    n_heavy_atoms  = 0
    for sub_SMILES in SMILES.split("."):
        n_heavy_atoms = max(n_heavy_atoms, fragmenter.get_heavy_atom_count(Chem.MolFromSmiles(sub_SMILES)))
    
    if n_heavy_atoms <= 20:
        right_size_for_complete_fragmenter2.append(inchikey)
        fragmentation, success = complete_fragmenter.fragment(SMILES)
        if success:
            complete_fragmenter_sorted_fragmented.append(inchikey)
            if is_fragmentation_equal_to_other_fragmentation(fragmentation, fragmentation_reference_DB):
                complete_fragmenter_sorted_fragmented_and_equal_to_reference_DB.append(inchikey)
        
        log_structure_results(f_complete, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB)
    else:
        log_structure_results(f_complete, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB, 'Structure was skipped because it is larger than 20 atoms.')

f_simple.close()
f_complete.close()
  
print ''
print 'N_structures(simple): ' + str(len(reference_DB))
print 'N_fragmented(simple): ' + str(len(simple_fragmenter_sorted_fragmented)) + "(" + str((1.0 * len(simple_fragmenter_sorted_fragmented)) / len(reference_DB)) + ")"
print 'N_fragmented_and_equal(simple): ' + str(len(simple_fragmenter_sorted_fragmented_and_equal_to_reference_DB)) + "(" + str((1.0 * len(simple_fragmenter_sorted_fragmented_and_equal_to_reference_DB)) / len(reference_DB)) + ")"
print ''
print 'N_structures(complete):' + str(len(right_size_for_complete_fragmenter2))
print 'N_fragmented(complete): ' + str(len(complete_fragmenter_sorted_fragmented)) + "(" + str((1.0 * len(complete_fragmenter_sorted_fragmented)) / len(right_size_for_complete_fragmenter2)) + ")"
print 'N_fragmented_and_equal(complete): ' + str(len(complete_fragmenter_sorted_fragmented_and_equal_to_reference_DB)) + "(" + str((1.0 * len(complete_fragmenter_sorted_fragmented_and_equal_to_reference_DB)) / len(right_size_for_complete_fragmenter2)) + ")"
print ''
print ''  


# second step: try to fragent all from the component    
structures_DB = []
with open('structures_DB.csv') as f:
    for line in f.readlines():
        structures_DB.append(CSV_to_info(line))
        
combined_fragmenter = fragmenter(fragmentation_scheme, 'combined', 20, function_to_choose_fragmentation, 1)
combined_fragmenter.fragmentation_scheme_order = sorted_group_numbers 
combined_fragmenter.n_max_fragmentations_to_find = 1    

combined_fragmenter_sorted_fragmented = []
right_size_for_combined_fragmenter  = []
print '####################################################################'
print 'Fragmenting the structures database with the patterns sorted (combined algorithm)'
i_structure = 0
f_combined = open('structures_DB_combined_fragmentation_with_pattern_sorting_results.log','w+')
for inchikey, SMILES, pubchem_id, empty_fragmentation in structures_DB:  
 
    i_structure = i_structure + 1
    if i_structure % 2000 == 0:
        print '{:2.1f} .'.format((100.0 * i_structure) / len(structures_DB)),
    
    fragmentation, success = combined_fragmenter.fragment(SMILES)
    if success:
        combined_fragmenter_sorted_fragmented.append(inchikey)
        
    n_heavy_atoms  = 0
    for sub_SMILES in SMILES.split("."):
        n_heavy_atoms = max(n_heavy_atoms, fragmenter.get_heavy_atom_count(Chem.MolFromSmiles(sub_SMILES)))
    
    if n_heavy_atoms <= 20:
        right_size_for_combined_fragmenter.append(inchikey)
        log_structure_results(f_combined, pubchem_id, SMILES, inchikey, success, fragmentation, {})
    else:
        if success:
            log_structure_results(f_combined, pubchem_id, SMILES, inchikey, success, fragmentation, {})
        else:
            log_structure_results(f_combined, pubchem_id, SMILES, inchikey, success, fragmentation, {}, 'Structure was skipped because it is larger than 20 atoms.')
          
  
f_combined.close()
    
print ''
print 'N_structures(simple): ' + str(len(right_size_for_combined_fragmenter))
print 'N_fragmented(simple): ' + str(len(combined_fragmenter_sorted_fragmented)) + "(" + str((1.0 * len(combined_fragmenter_sorted_fragmented)) / len(right_size_for_combined_fragmenter)) + ")"
print ''
print '####################################################################'