# -*- coding: utf-8 -*-
''' Class for fragmenting molecules into molecular subgroups

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

    
    
class fragmenter:
    # tested with Python 3.6.8 and RDKit version 2017.09.3
    
    # import dependencies
    try:
        import rdkit as __rdkit
        from rdkit import Chem as __Chem
        import marshal as __marshal
        import regex as __regex

    except:
        raise Exception('rdkit, marshal and regex have to be installed.')
        
        
    if __rdkit.rdBase.rdkitVersion != '2017.09.3':
        raise Exception('In this code, the SMARTS are \
analyzed for their properties. Unfortunately different rdkit versions give different values. \
This code has only been tested on version \'2017.09.3\'.')
    # for more details have a look at the comments in the function get_mol_with_properties_from_SMARTS
    
                    
    # get a molecule from a SMARTS with the properties necessary to calculate
    # more complex properties
    @classmethod
    def get_mol_from_SMARTS(cls, smarts):
        mol = cls.__Chem.MolFromSmarts(smarts)
        mol.UpdatePropertyCache(strict = False)
        cls.__Chem.GetSymmSSSR(mol)
        return mol   
    
    # get a molecule from a SMARTS and calculate complex properties
    # this function could be improved by someone who knows RDkit better than me
	# 
	# Unfortunately it only works corrrectly with version 2017.09.3 of RDKit,
	# for more details see the following issues:
	# https://github.com/rdkit/rdkit/issues/2448
	# https://github.com/rdkit/rdkit/issues/1978
	# 
	# if anyone finds a better way of calculating the number of
	# bonds a SMARTS fragment has, please contact me or send a merge request
    #
    # if you are using conda you may install it in a new environemnt using the following command: conda install -c rdkit rdkit=2017.09.3
	
    @classmethod
    def get_mol_with_properties_from_SMARTS(cls, SMARTS):
        if cls.__rdkit.rdBase.rdkitVersion != '2017.09.3':
            print('#################WARNING########################')
            print('In this code, the SMARTS are \
analyzed programmatically for their properties. Unfortunately different rdkit versions give different values. \
This code has been developed with version \'2017.09.3\'.')
            # if you are using conda you may install it in a new environemnt using the following command: conda install -c rdkit rdkit=2017.09.3
            print('################################################')
                  
        mol_SMARTS = fragmenter.get_mol_from_SMARTS(SMARTS)
        
        # some cleaning up
        # conditional SMARTS are not supported
        SMARTS = SMARTS.replace('-,:', '')
        conditional_statement_exists = SMARTS.count(',') > 0
        if conditional_statement_exists:
            print (SMARTS)
            raise ValueError('Algorithm can\'t handle conditional SMARTS. Please use a list of SMARTS in the fragmentation scheme to mimic conditional SMARTS.')
        
        n_atoms_defining_SMARTS = 0
        n_available_bonds = 0
        n_atoms_with_available_bonds = 0
        n_hydrogens_this_molecule = 0
        n_atoms_in_ring = 0
        n_carbons = 0
        found_atoms = []
        
        # iterate over atoms to get their properties
        for atom in mol_SMARTS.GetAtoms():
            SMARTS_atom = atom.GetSmarts()
            n_atoms_defining_SMARTS += 1
            n_available_bonds_this_atom = 0
            
            if atom.GetSymbol() == '*':
                matches = cls.__regex.finditer("\$?\(([^()]|(?R))*\)", SMARTS_atom)
                # if it is a recursive SMARTS
                if matches is not None:
                    for match in matches:
                        m = match.group(0)
                        
                        found_atoms.append(m[2:-1])
                        mol_SMARTS2 = fragmenter.get_mol_with_properties_from_SMARTS(m[2:-1])
                        n_atoms_defining_SMARTS += mol_SMARTS2.GetUnsignedProp('n_atoms_defining_SMARTS')
                        
                        first_atom = mol_SMARTS2.GetAtomWithIdx(0)
                        n_hydrogens = first_atom.GetUnsignedProp('n_hydrogens')
                        
                        n_available_bonds_this_atom = first_atom.GetTotalValence() - n_hydrogens
                        
                        atom.SetUnsignedProp('n_hydrogens', n_hydrogens)
                        atom.SetUnsignedProp('n_available_bonds', n_available_bonds_this_atom)
                        
                        if first_atom.GetAtomicNum() == 6:
                            n_carbons += 1
                            
                    if len(found_atoms) == 1:
                        n_atoms_defining_SMARTS -= 1
                    elif len(found_atoms) > 1:
                        raise ValueError('Algorithm can\'t handle SMARTS with 2 levels of recursion')
                    else:
                        atom.SetUnsignedProp('n_hydrogens', 0)
                        atom.SetUnsignedProp('n_available_bonds', 0)
            else:
                
                # get number of hydrogens from SMARTS
                n_available_bonds_this_atom = atom.GetImplicitValence()
                
                n_hydrogens = 0
                match = cls.__regex.findall('AtomHCount\s+(\d+)\s*', atom.DescribeQuery())
                if match:
                    n_hydrogens = int(match[0])
                    
                n_available_bonds_this_atom -= n_hydrogens
                    
                atom.SetUnsignedProp('n_hydrogens', n_hydrogens)
                
                if n_available_bonds_this_atom < 0:
                    n_available_bonds_this_atom = 0
                
                atom.SetUnsignedProp('n_available_bonds', n_available_bonds_this_atom)
                
                if atom.GetAtomicNum() == 6:
                    n_carbons += 1
                    
                is_within_ring = False
                match = cls.__regex.findall('AtomInNRings\s+(-?\d+)\s+(!?=)\s*', atom.DescribeQuery())
                if match:
                     match = match[0]
                     n_rings = int(match[0])
                     comparison_sign = match[1]
                     
                     if n_rings == -1: # if number of rings was not defined 
                         is_within_ring = comparison_sign == '='
                     else: # if number of rings was defined
                         if comparison_sign == '=':
                             is_within_ring  = n_rings != 0
                         else:
                             is_within_ring  = n_rings == 0
                        
                if atom.GetIsAromatic() or is_within_ring:
                    n_atoms_in_ring += 1
                    
            n_available_bonds +=n_available_bonds_this_atom
            
            n_hydrogens_this_molecule += atom.GetUnsignedProp('n_hydrogens')
            
            if n_available_bonds_this_atom > 0:
                n_atoms_with_available_bonds += 1

        
        # find whether the SMARTS is simple
        atom_with_valence_one_on_carbon = cls.__Chem.MolFromSmarts('[*;v1][#6]')
        atom_with_valence_one_on_carbon.UpdatePropertyCache()
        
        atom_with_valence_one_on_excluding_carbon = cls.__Chem.MolFromSmarts('[$([*;v1][#6])]')
        atom_with_valence_one_on_excluding_carbon.UpdatePropertyCache()

        is_simple_atom_on_c = n_atoms_with_available_bonds == 1 and \
                                    (mol_SMARTS.HasSubstructMatch(atom_with_valence_one_on_carbon) or \
                                     mol_SMARTS.HasSubstructMatch(atom_with_valence_one_on_excluding_carbon))
        
        atom_with_valence_one = cls.__Chem.MolFromSmarts('[*;v1;!#1]')
        atom_with_valence_one.UpdatePropertyCache()
        is_simple_atom = n_atoms_with_available_bonds == 1 and (mol_SMARTS.HasSubstructMatch(atom_with_valence_one))

        
        if len(found_atoms) > 0:
            if len(found_atoms) == 1:
                sub_mol_SMARTS = cls.__Chem.MolFromSmarts(found_atoms[0])
                sub_mol_SMARTS.UpdatePropertyCache()
                is_simple_atom_on_c = (sub_mol_SMARTS.HasSubstructMatch(atom_with_valence_one_on_carbon) or \
                                       sub_mol_SMARTS.HasSubstructMatch(atom_with_valence_one_on_excluding_carbon))     
                is_simple_atom = n_atoms_defining_SMARTS == 1 and sub_mol_SMARTS.HasSubstructMatch(atom_with_valence_one)
            elif len(found_atoms) > 1:
                raise ValueError('Algorithm can\'t handle SMARTS with 2 recursive SMARTS')
          
        # set the gathered properties
        mol_SMARTS.SetUnsignedProp('n_hydrogens', n_hydrogens_this_molecule)
        mol_SMARTS.SetUnsignedProp('n_carbons', n_carbons)
        mol_SMARTS.SetUnsignedProp('n_available_bonds', n_available_bonds)
        mol_SMARTS.SetUnsignedProp('n_atoms_defining_SMARTS', n_atoms_defining_SMARTS)                
        mol_SMARTS.SetUnsignedProp('n_atoms_with_available_bonds', n_atoms_with_available_bonds)
        mol_SMARTS.SetBoolProp('has_atoms_in_ring', n_atoms_in_ring > 0)
        mol_SMARTS.SetBoolProp('is_simple_atom', is_simple_atom)
        mol_SMARTS.SetBoolProp('is_simple_atom_on_c', is_simple_atom_on_c)
        mol_SMARTS.SetUnsignedProp('n_double_bonds', len(mol_SMARTS.GetSubstructMatches(cls.__Chem.MolFromSmarts("*=*"))))
        mol_SMARTS.SetUnsignedProp('n_triple_bonds', len(mol_SMARTS.GetSubstructMatches(cls.__Chem.MolFromSmarts("*#*"))))
        
        return mol_SMARTS      

    # this function does a substructure match and then checks whether the match 
    # is adjacent to previous matches and/or checks if the hydrogen number is correct
    @classmethod
    def get_substruct_matches(cls, mol_searched_for, mol_searched_in, atomIdxs_to_which_new_matches_have_to_be_adjacent, check_for_hydrogens = False):
        
        valid_matches = []
        
        if mol_searched_in.GetNumAtoms() >= mol_searched_for.GetNumAtoms():
            matches = mol_searched_in.GetSubstructMatches(mol_searched_for)
            
            if matches:
                for match in matches:
                    all_hydrogens_OK = True
                    
                    if check_for_hydrogens:
                        # following lines are a workaround for the fact that SMARTS
                        # matching SMARTS is not working completely correctly as the number
                        # of hydrogens is ignored in some cases by RDkit
                        # 
                    	# for more details see the following issues:
                    	# https://github.com/rdkit/rdkit/issues/2448
                    	# https://github.com/rdkit/rdkit/issues/1978
                        
                        for i in range(mol_searched_for.GetNumAtoms()):
                            atom_mol_searched_for = mol_searched_for.GetAtomWithIdx(i)
                            atom_mol_searched_in = mol_searched_in.GetAtomWithIdx(match[i])
                            
                            # if mol_searched_in is SMARTS
                            if atom_mol_searched_in.HasProp('n_hydrogens'):
                                if atom_mol_searched_for.GetUnsignedProp('n_hydrogens') > atom_mol_searched_in.GetUnsignedProp('n_hydrogens'):
                                    all_hydrogens_OK = False 
                                    break
                            # if mol_searched_in is SMILES
                            else:
                                break
                        
                    if all_hydrogens_OK:
                        add_this_match = True
                        if len(atomIdxs_to_which_new_matches_have_to_be_adjacent) > 0:
                            add_this_match = False
                            
                            for i in match:
                                for neighbor in mol_searched_in.GetAtomWithIdx(i).GetNeighbors():
                                    if neighbor.GetIdx() in atomIdxs_to_which_new_matches_have_to_be_adjacent:
                                        add_this_match = True
                                        break
                                
                        if add_this_match:
                            valid_matches.append(match)
                
        return valid_matches
    
    # this dunction is to avoid counting heavier versions of hydrogen as heavy atom
    @classmethod
    def get_heavy_atom_count(cls, mol):
        heavy_atom_count = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 1:
                heavy_atom_count += 1
        
        return heavy_atom_count
    
    def __init__(self, fragmentation_scheme = {}, algorithm = '', n_atoms_cuttoff = -1, function_to_choose_fragmentation = False, n_max_fragmentations_to_find = -1):
    
        if not type(fragmentation_scheme) is dict:
            raise TypeError('fragmentation_scheme must be a dctionary with integers as keys and either strings or list of strings as values.')
            
        if len(fragmentation_scheme) == 0:
            raise ValueError('fragmentation_scheme must be provided.')
        
        if not algorithm in ['simple', 'complete', 'combined']:
            raise ValueError('Algorithm must be either simple ,complete or combined.')
            
        if algorithm == 'simple':
            if n_max_fragmentations_to_find != -1:
                raise ValueError('Setting n_max_fragmentations_to_find only makes sense with complete or combined algorithm.')
        
        self.algorithm = algorithm
        
        if algorithm in ['complete', 'combined']:
            if n_atoms_cuttoff == -1:
                raise ValueError('n_atoms_cuttoff needs to be specified for complete or combined algorithms.')
                
            if function_to_choose_fragmentation == False:
                raise ValueError('function_to_choose_fragmentation needs to be specified for complete or combined algorithms.')
                
            if not callable(function_to_choose_fragmentation):
                raise TypeError('function_to_choose_fragmentation needs to be a function.')
                
            if n_max_fragmentations_to_find != -1:
                if n_max_fragmentations_to_find < 1:
                    raise ValueError('n_max_fragmentations_to_find has to be 1 or higher.')
            
        self.n_max_fragmentations_to_find = n_max_fragmentations_to_find
        
        self.n_atoms_cuttoff = n_atoms_cuttoff
        
        self.fragmentation_scheme = fragmentation_scheme
        
        self.function_to_choose_fragmentation = function_to_choose_fragmentation
        
        # create a lookup dictionaries to faster finding a group number and the 
        # respective pattern with its properties for a specific SMARTS
        self._fragmentation_scheme_group_number_lookup = {}
        self._fragmentation_scheme_pattern_lookup = {}
        self.fragmentation_scheme_order = []
        for group_number, list_SMARTS in fragmentation_scheme.items():
            
            self.fragmentation_scheme_order.append(group_number)
            
            if type(list_SMARTS) is not list:
                list_SMARTS = [list_SMARTS]
                
            for SMARTS in list_SMARTS:
                if SMARTS != '':
                    self._fragmentation_scheme_group_number_lookup[SMARTS] = group_number
                    
                    mol_SMARTS = fragmenter.get_mol_with_properties_from_SMARTS(SMARTS)
                    self._fragmentation_scheme_pattern_lookup[SMARTS] = mol_SMARTS
        
        

        # create a lookup dictionaries to faster finding of parent pattern
        # for a specific SMARTS
        self._parent_pattern_lookup = {}
        for SMARTS1, mol_SMARTS1 in self._fragmentation_scheme_pattern_lookup.items():
            parent_patterns_of_SMARTS1 = []
            
            for SMARTS2, mol_SMARTS2 in self._fragmentation_scheme_pattern_lookup.items():
                if SMARTS1 != SMARTS2:
                    if mol_SMARTS2.GetNumAtoms() >= mol_SMARTS1.GetNumAtoms():
                        matches = fragmenter.get_substruct_matches(mol_SMARTS1, mol_SMARTS2, set(), True)
                        
                        if matches:
                            parent_patterns_of_SMARTS1.append(SMARTS2)
            
            mol_SMARTS1.SetBoolProp('has_parent_pattern', len(parent_patterns_of_SMARTS1) > 0)

            self._parent_pattern_lookup[SMARTS1] = parent_patterns_of_SMARTS1


    def fragment(self, SMILES):
        
        is_valid_SMILES = fragmenter.__Chem.MolFromSmiles(SMILES) is not None
        
        if not is_valid_SMILES:
            raise ValueError('Following SMILES is not valid: ' + SMILES)
        
        # handle mixtures
        if SMILES.count('.') > 0:
            list_SMILES = SMILES.split('.')
        else:
            list_SMILES = [SMILES]
            
        # iterate over all separated SMILES
        success = False
        fragmentation = {}
        for SMILES in list_SMILES:
            temp_fragmentation, success = self.__get_fragmentation(SMILES)
    
            for SMARTS, matches in temp_fragmentation.items():
                group_number = self._fragmentation_scheme_group_number_lookup[SMARTS]
                
                if not group_number in fragmentation:
                    fragmentation[group_number] = 0
                
                fragmentation[group_number] += len(matches)   
                
            if not success:
                break
        
        return fragmentation, success
        
    def __get_fragmentation(self, SMILES):
        
        success = False
        fragmentation = {}
        # use simple fragmentation algorithm
        if self.algorithm in ['simple', 'combined']:
            fragmentation, success = self.__simple_fragmentation(SMILES)
        
        if success:
            return fragmentation, success
        
        if self.algorithm in ['complete', 'combined']:
            fragmentations, success = self.__complete_fragmentation(SMILES)
            
            if success:
                fragmentation = self.function_to_choose_fragmentation(fragmentations)
        
        return fragmentation, success
    
    def __simple_fragmentation(self, SMILES):
        mol_SMILES = self.__Chem.MolFromSmiles(SMILES)
    
        heavy_atom_count = fragmenter.get_heavy_atom_count(mol_SMILES)
    
        success = False
        fragmentation = {}
        
        fragmentation, atomIdxs_included_in_fragmentation = self.__search_non_overlapping_solution(mol_SMILES, {}, set(), set())
        success = len(atomIdxs_included_in_fragmentation) == heavy_atom_count
        
        # if not success, clean up molecule and search again
        level = 1
        while not success:
            fragmentation_so_far , atomIdxs_included_in_fragmentation_so_far = fragmenter.__clean_molecule_surrounding_unmatched_atoms(mol_SMILES, fragmentation, atomIdxs_included_in_fragmentation, level)
            level += 1
            
            if len(atomIdxs_included_in_fragmentation_so_far) == 0:
                break
            
            fragmentation_so_far, atomIdxs_included_in_fragmentation_so_far = self.__search_non_overlapping_solution(mol_SMILES, fragmentation_so_far, atomIdxs_included_in_fragmentation_so_far, atomIdxs_included_in_fragmentation_so_far)
            
            success = len(atomIdxs_included_in_fragmentation_so_far) == heavy_atom_count
            
            if success:
                fragmentation = fragmentation_so_far
            
        return fragmentation, success
    
    
    def __search_non_overlapping_solution(self, mol_searched_in, fragmentation, atomIdxs_included_in_fragmentation, atomIdxs_to_which_new_matches_have_to_be_adjacent, search_for_parent_patterns = True):
        
        n_atomIdxs_included_in_fragmentation = len(atomIdxs_included_in_fragmentation) - 1
        
        while n_atomIdxs_included_in_fragmentation != len(atomIdxs_included_in_fragmentation):
            n_atomIdxs_included_in_fragmentation = len(atomIdxs_included_in_fragmentation)
            
             
            for group_number in self.fragmentation_scheme_order:
                list_SMARTS = self.fragmentation_scheme[group_number]
                if type(list_SMARTS) is not list:
                    list_SMARTS = [list_SMARTS]
                
                for SMARTS in list_SMARTS:
                    if SMARTS != "":  
                        fragmentation, atomIdxs_included_in_fragmentation = self.__get_next_non_overlapping_match(mol_searched_in, SMARTS, fragmentation, atomIdxs_included_in_fragmentation, atomIdxs_to_which_new_matches_have_to_be_adjacent, search_for_parent_patterns)

        return fragmentation, atomIdxs_included_in_fragmentation

    def __get_next_non_overlapping_match(self, mol_searched_in, SMARTS, fragmentation, atomIdxs_included_in_fragmentation, atomIdxs_to_which_new_matches_have_to_be_adjacent, search_for_parent_patterns):
        
        mol_searched_for = self._fragmentation_scheme_pattern_lookup[SMARTS]
        
        if search_for_parent_patterns:
            if mol_searched_for.GetBoolProp('has_parent_pattern'):
                for parent_SMARTS in self._parent_pattern_lookup[SMARTS]:
                    fragmentation, atomIdxs_included_in_fragmentation = self.__get_next_non_overlapping_match(mol_searched_in, parent_SMARTS, fragmentation, atomIdxs_included_in_fragmentation, atomIdxs_to_which_new_matches_have_to_be_adjacent, search_for_parent_patterns)
       
        if atomIdxs_to_which_new_matches_have_to_be_adjacent:
            matches = fragmenter.get_substruct_matches(mol_searched_for, mol_searched_in, atomIdxs_to_which_new_matches_have_to_be_adjacent)
        else:
            matches = fragmenter.get_substruct_matches(mol_searched_for, mol_searched_in, set())
        
        if matches:
            for match in matches:
                all_atoms_of_new_match_are_unassigned = atomIdxs_included_in_fragmentation.isdisjoint(match)
        
                if all_atoms_of_new_match_are_unassigned:
                    if not SMARTS in fragmentation:
                        fragmentation[SMARTS] = []
                        
                    fragmentation[SMARTS].append(match)
                    atomIdxs_included_in_fragmentation.update(match)   
                
        return fragmentation, atomIdxs_included_in_fragmentation

    @classmethod
    def __clean_molecule_surrounding_unmatched_atoms(cls, mol_searched_in, fragmentation, atomIdxs_included_in_fragmentation, level):
    
        for i in range(0, level):
            
            atoms_missing = set(range(0, fragmenter.get_heavy_atom_count(mol_searched_in))).difference(atomIdxs_included_in_fragmentation)
                        
            new_fragmentation = fragmenter.__marshal.loads(fragmenter.__marshal.dumps(fragmentation))
            
            for atomIdx in atoms_missing:
                for neighbor in mol_searched_in.GetAtomWithIdx(atomIdx).GetNeighbors():
                    for smart, atoms_found in fragmentation.items():
                        for atoms in atoms_found:
                            if neighbor.GetIdx() in atoms:
                                if smart in new_fragmentation:
                                    if new_fragmentation[smart].count(atoms) > 0:
                                        new_fragmentation[smart].remove(atoms)
                                
                        if smart in new_fragmentation:
                            if len(new_fragmentation[smart]) == 0:
                                new_fragmentation.pop(smart)
                                
                                
            new_atomIdxs_included_in_fragmentation = set()
            for i in new_fragmentation.values():
                for j in i:
                    new_atomIdxs_included_in_fragmentation.update(j)
                    
            atomIdxs_included_in_fragmentation = new_atomIdxs_included_in_fragmentation
            fragmentation = new_fragmentation
            
        return fragmentation, atomIdxs_included_in_fragmentation

       
    def __complete_fragmentation(self, SMILES):
        mol_SMILES = self.__Chem.MolFromSmiles(SMILES)
    
        heavy_atom_count = fragmenter.get_heavy_atom_count(mol_SMILES)
        
        if heavy_atom_count > self.n_atoms_cuttoff:
            return {}, False
        
        completed_fragmentations = []
        groups_leading_to_incomplete_fragmentations = []        
        completed_fragmentations, groups_leading_to_incomplete_fragmentations, incomplete_fragmentation_found = self.__get_next_non_overlapping_adjacent_match_recursively(mol_SMILES, heavy_atom_count, completed_fragmentations, groups_leading_to_incomplete_fragmentations, {}, set(), set(), self.n_max_fragmentations_to_find)
        success = len(completed_fragmentations) > 0
        
        return completed_fragmentations, success
        
    def __get_next_non_overlapping_adjacent_match_recursively(self, mol_searched_in, heavy_atom_count, completed_fragmentations, groups_leading_to_incomplete_fragmentations, fragmentation_so_far, atomIdxs_included_in_fragmentation_so_far, atomIdxs_to_which_new_matches_have_to_be_adjacent, n_max_fragmentations_to_find = -1):
      
        n_completed_fragmentations = len(completed_fragmentations)
        incomplete_fragmentation_found = False
        complete_fragmentation_found = False
        
        if len(completed_fragmentations) == n_max_fragmentations_to_find:
            return completed_fragmentations, groups_leading_to_incomplete_fragmentations, incomplete_fragmentation_found
                    
                
        for group_number in self.fragmentation_scheme_order:
            list_SMARTS = self.fragmentation_scheme[group_number]
            
            if complete_fragmentation_found:
                break
            
            if type(list_SMARTS) is not list:
                list_SMARTS = [list_SMARTS]
            
            for SMARTS in list_SMARTS:
                if complete_fragmentation_found:
                    break
                
                if SMARTS != "":  
                    matches = fragmenter.get_substruct_matches(self._fragmentation_scheme_pattern_lookup[SMARTS], mol_searched_in, atomIdxs_included_in_fragmentation_so_far)
                    
                    for match in matches:
                        
                        # only allow non-overlapping matches
                        all_atoms_are_unassigned = atomIdxs_included_in_fragmentation_so_far.isdisjoint(match)
                        if not all_atoms_are_unassigned:
                            continue
                        
                        # only allow matches that do not contain groups leading to incomplete matches
                        for groups_leading_to_incomplete_fragmentation in groups_leading_to_incomplete_fragmentations:
                            if fragmenter.__is_fragmentation_subset_of_other_fragmentation(groups_leading_to_incomplete_fragmentation, fragmentation_so_far):
                                return completed_fragmentations, groups_leading_to_incomplete_fragmentations, incomplete_fragmentation_found
                        
                        # only allow matches that will lead to new fragmentations
                        use_this_match = True
                        n_found_groups = len(fragmentation_so_far)
                        
                        for completed_fragmentation in completed_fragmentations:
                            
                            if not SMARTS in completed_fragmentation:
                                continue
                            
                            if n_found_groups == 0:
                                use_this_match = not fragmenter.__is_match_contained_in_fragmentation(match, SMARTS, completed_fragmentation)
                            else:
                                if fragmenter.__is_fragmentation_subset_of_other_fragmentation(fragmentation_so_far, completed_fragmentation):
                                    use_this_match = not fragmenter.__is_match_contained_in_fragmentation(match, SMARTS, completed_fragmentation)
                            
                            if not use_this_match:
                                break
                                
                        if not use_this_match:
                            continue
                        
                        # make a deepcopy here, otherwise the variables are modified down the road
                        # marshal is used here because it works faster than copy.deepcopy
                        this_SMARTS_fragmentation_so_far = fragmenter.__marshal.loads(fragmenter.__marshal.dumps(fragmentation_so_far))
                        this_SMARTS_atomIdxs_included_in_fragmentation_so_far = atomIdxs_included_in_fragmentation_so_far.copy()
                        
                        if not SMARTS in this_SMARTS_fragmentation_so_far:
                            this_SMARTS_fragmentation_so_far[SMARTS] = []
                            
                        this_SMARTS_fragmentation_so_far[SMARTS].append(match)
                        this_SMARTS_atomIdxs_included_in_fragmentation_so_far.update(match)
                        
                        # only allow matches that do not contain groups leading to incomplete matches
                        for groups_leading_to_incomplete_match in groups_leading_to_incomplete_fragmentations:
                            if fragmenter.__is_fragmentation_subset_of_other_fragmentation(groups_leading_to_incomplete_match, this_SMARTS_fragmentation_so_far):
                                use_this_match = False
                                break
                            
                        if not use_this_match:
                            continue
                        
                        # if the complete molecule has not been fragmented, continue to do so
                        if len(this_SMARTS_atomIdxs_included_in_fragmentation_so_far) < heavy_atom_count:
                            completed_fragmentations, groups_leading_to_incomplete_fragmentations, incomplete_fragmentation_found = self.__get_next_non_overlapping_adjacent_match_recursively(mol_searched_in, heavy_atom_count, completed_fragmentations, groups_leading_to_incomplete_fragmentations, this_SMARTS_fragmentation_so_far, this_SMARTS_atomIdxs_included_in_fragmentation_so_far, this_SMARTS_atomIdxs_included_in_fragmentation_so_far, n_max_fragmentations_to_find)
                            break
                        
                        # if the complete molecule has been fragmented, save and return
                        if len(this_SMARTS_atomIdxs_included_in_fragmentation_so_far) == heavy_atom_count:                 
                            completed_fragmentations.append(this_SMARTS_fragmentation_so_far)
                            complete_fragmentation_found = True
                            break
                        
        # if until here no new fragmentation was found check whether an incomplete fragmentation was found
        if n_completed_fragmentations == len(completed_fragmentations):      
            
            if not incomplete_fragmentation_found:
                
                incomplete_matched_groups = {}
                
                if len(atomIdxs_included_in_fragmentation_so_far) > 0:
                    unassignes_atom_idx = set(range(0, heavy_atom_count)).difference(atomIdxs_included_in_fragmentation_so_far)
                    for atom_idx in unassignes_atom_idx:
                        neighbor_atoms_idx = [i.GetIdx() for i in mol_searched_in.GetAtomWithIdx(atom_idx).GetNeighbors()]
                        
                        for neighbor_atom_idx in neighbor_atoms_idx:
                            for found_smarts, found_matches in fragmentation_so_far.items():
                                for found_match in found_matches:
                                    if neighbor_atom_idx in found_match:
                                        if not found_smarts in incomplete_matched_groups:
                                            incomplete_matched_groups[found_smarts] = []
                                            
                                        if found_match not in incomplete_matched_groups[found_smarts]:
                                            incomplete_matched_groups[found_smarts].append(found_match)
                                    
                    is_subset_of_groups_already_found = False
                    indexes_to_remove = []
                    ind = 0
                    
                    # remove groups that are parent to the currently found groups
                    for groups_leading_to_incomplete_match in groups_leading_to_incomplete_fragmentations:
                        is_subset_of_groups_already_found = fragmenter.__is_fragmentation_subset_of_other_fragmentation(incomplete_matched_groups, groups_leading_to_incomplete_match)
                        if is_subset_of_groups_already_found:
                            indexes_to_remove.append(ind)
                            
                        ind += 1
                        
                    for index in sorted(indexes_to_remove, reverse=True):
                        del groups_leading_to_incomplete_fragmentations[index]
                        
                    groups_leading_to_incomplete_fragmentations.append(incomplete_matched_groups)
                    groups_leading_to_incomplete_fragmentations = sorted(groups_leading_to_incomplete_fragmentations, key = len)
                    
                    incomplete_fragmentation_found =  True
    
        return completed_fragmentations, groups_leading_to_incomplete_fragmentations, incomplete_fragmentation_found
    
    @classmethod
    def __is_fragmentation_subset_of_other_fragmentation(cls, fragmentation, other_fragmentation):
        n_found_groups = len(fragmentation)
        n_found_other_groups = len(other_fragmentation)
        
        if n_found_groups == 0:
            return False
            
        if n_found_other_groups < n_found_groups:
            return False
        
        n_found_SMARTS_that_are_subset = 0
        for found_SMARTS, found_matches in fragmentation.items():
            if found_SMARTS in other_fragmentation:
                found_matches_set = set(frozenset(i) for i in fragmentation[found_SMARTS])
                found_other_matches_set =  set(frozenset(i) for i in other_fragmentation[found_SMARTS])
                
                if found_matches_set.issubset(found_other_matches_set):
                    n_found_SMARTS_that_are_subset += 1
            else:
                return False
            
        return n_found_SMARTS_that_are_subset == n_found_groups
    
    @classmethod
    def __is_match_contained_in_fragmentation(cls, match, SMARTS, fragmentation):
        if not SMARTS in fragmentation:
            return False
            
        found_matches_set = set(frozenset(i) for i in fragmentation[SMARTS])
        match_set = set(match)
        
        return match_set in found_matches_set