# -*- coding: gbk -*-
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem
from rdkit.Chem import rdmolops
from rdkit.Chem import rdDepictor
from rdkit.Chem import EditableMol
import math
from collections import deque
from random import shuffle
import re

def process_molecule(mol):
    smiles = Chem.MolToSmiles(mol)
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    smiles = Chem.MolToSmiles(mol)
    return mol, smiles

def get_atom_groups(mol):
    # ��ȡ���������л�����Ϣ
    ssr = Chem.GetSymmSSSR(mol)

    # �ϲ���һ������ԭ���غϵĻ�
    atom_groups = []
    merged_rings = []
    for ring in ssr:
        merged = False
        for i, merged_ring in enumerate(merged_rings):
            if len(set(ring) & set(merged_ring)) >= 1:
                merged_rings[i] = tuple(set(ring) | set(merged_ring))
                merged = True
                break
        if not merged:
            merged_rings.append(ring)
            
    while True:
        merged = False
        for i, ring1 in enumerate(merged_rings):
            for j, ring2 in enumerate(merged_rings[i+1:]):
                if len(set(ring1) & set(ring2)) >= 1:
                    merged_rings[i] = tuple(set(ring1) | set(ring2))
                    merged_rings.pop(i+j+1)
                    merged = True
                    break
            if merged:
                break
        if not merged:
            break
        
        
    # ����ÿ����
    for ring_idx, ring in enumerate(merged_rings):
        # ����ÿ�����е�ÿ��ԭ��
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # ���ԭ���Ƿ��뵥��ԭ����˫������������
            for bond in atom.GetBonds():
                if bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
                    other_atom_idx = bond.GetOtherAtomIdx(atom_idx)
                    # ����ԭ���Ƿ�����һ������
                    in_other_ring = False
                    for other_ring_idx, other_ring in enumerate(merged_rings):
                        if other_ring != ring and other_atom_idx in other_ring:
                            # ���������ϲ�Ϊһ���µĻ�
                            ring = tuple(set(ring) | set(other_ring))
                            merged_rings[ring_idx] = ring
                            del merged_rings[other_ring_idx]
                            in_other_ring = True
                            break
                    if not in_other_ring:
                        # ����ԭ����ӵ�����
                        ring = tuple(set(ring) | {other_atom_idx})
                        merged_rings[ring_idx] = ring

    atom_groups = [list(ring) for ring in merged_rings]


    non_ring_atoms = set(range(mol.GetNumAtoms())) - set().union(*merged_rings)
    non_ring_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in non_ring_atoms and bond.GetEndAtomIdx() in non_ring_atoms:
            #if bond.GetBondType() != Chem.BondType.SINGLE
            if bond.GetBondType() != Chem.BondType.SINGLE or bond.GetBeginAtom().GetSymbol() != 'C' or bond.GetEndAtom().GetSymbol() != 'C' :
                non_ring_bonds.append(bond.GetIdx())
    
    for atom_idx in non_ring_atoms:
        atom_groups.append(set([atom_idx]))
    
    for bond_idx in non_ring_bonds:
        bond = mol.GetBondWithIdx(bond_idx)
        begin_atom_idx = bond.GetBeginAtomIdx()
        end_atom_idx = bond.GetEndAtomIdx()
        begin_group = None
        end_group = None
        for group in atom_groups:
            if begin_atom_idx in group:
                begin_group = group
            if end_atom_idx in group:
                end_group = group
        if begin_group is not None and end_group is not None and begin_group != end_group:
            new_group = begin_group.union(end_group)
            atom_groups.remove(begin_group)
            atom_groups.remove(end_group)
            atom_groups.append(new_group)

    return atom_groups


def mark_atoms(mol, group):
    # ���������е�ԭ��
    for atom_idx in group:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in group:
                # Ϊԭ������Զ�����1
                neighbor.SetAtomMapNum(1)

    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == 1:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in group:
                    # ����ھ�ԭ���ڻ����У���Ϊ������Զ�����2
                    neighbor.SetAtomMapNum(2)
                       
def get_groups_mol(mol, left_atoms):
    left_groups = set(range(mol.GetNumAtoms())) - set().union(left_atoms)
    left_groups_bonds = [bond.GetIdx() for bond in mol.GetBonds() if bond.GetBeginAtomIdx() in left_groups and bond.GetEndAtomIdx() in left_groups]
    left_groups_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(left_groups), bondsToUse=left_groups_bonds)
    groups_mol = Chem.MolFromSmiles(left_groups_smiles, sanitize=False)
    return groups_mol

def bfs(mol, start_atom):
    # ���������������
    visited = set()
    order = []
    queue = deque([(start_atom, 0)])
    while queue:
        atom, level = queue.popleft()
        if atom.GetIdx() not in visited:
            visited.add(atom.GetIdx())
            if len(order) <= level:
                order.append([])
            order[level].append((atom.GetSymbol(), get_bond_types(atom)))
            #order[level].append((atom.GetSymbol()))
            neighbors = list(atom.GetNeighbors())
            '''
            all_neighbors_visited = True
            for neighbor in neighbors:
                if neighbor.GetIdx() not in visited:
                    all_neighbors_visited = False
                    break
            if all_neighbors_visited:
                order[level].append(("0","0"))
            '''
            for neighbor in neighbors:
                queue.append((neighbor, level + 1))
    for level in order:
        level.sort()
    return order
    
def get_bond_types(atom):
    bond_types = []
    for bond in atom.GetBonds():
        bond_types.append(str(bond.GetBondType()))
    bond_types.sort()
    return bond_types
    
def get_atom_info(atom: rdchem.Atom) -> tuple:
    symbol = atom.GetSymbol()
    charge = atom.GetFormalCharge()
    hybridization = str(atom.GetHybridization())
    spin_multiplicity = atom.GetNumRadicalElectrons() + 1
    valence = atom.GetDegree()
    isotope = atom.GetIsotope()
    aromatic = atom.GetIsAromatic()

    # ��ȡԭ�����ڷ���
    mol = atom.GetOwningMol()
    
    ssr = Chem.GetSymmSSSR(mol)
    atom_groups = []
    merged_rings = []
    for ring in ssr:
        merged_rings.append(ring)
        
    ring_fragments_smiles = []
    atom_groups = [list(ring) for ring in merged_rings]
    for groups in atom_groups:
        if atom.GetIdx() in groups:
            smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(groups))
            test_mol = Chem.MolFromSmiles(smiles, sanitize=False)
            smiles = Chem.MolToSmiles(test_mol)
            ring_fragments_smiles.append(smiles)
    ring_fragments_smiles.sort()
    
    order = bfs(mol, atom)
            
    return (symbol, charge, hybridization, spin_multiplicity, valence, isotope, aromatic,ring_fragments_smiles,order)
    #return (symbol, charge, hybridization, spin_multiplicity, valence, mass, isotope, aromatic, tuple(neighbors),bond_types, bond_angles,ring_fragments_smiles)
   
def remove_groups(mol, atom_groups):
    iteration = 0
    successful_groups = []
    start_mol = None
    while True:
        iteration += 1
        removed = False
        for group in atom_groups:
            #������
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
            # ��¼����λ��
            connected_atoms = []
            connected_group_atoms = [] 
            # ��¼���ɾ���Ļ��ŷ��Ӷ���
            start_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(group))       
            start_mol = Chem.MolFromSmiles(start_smiles, sanitize=False)
            #���ϱ��
            mark_atoms(mol, group)

            left_atoms = set(range(mol.GetNumAtoms())) - set().union(group)
            if not left_atoms:
                return successful_groups, start_mol

            left_bonds = [bond.GetIdx() for bond in mol.GetBonds() if bond.GetBeginAtomIdx() in left_atoms and bond.GetEndAtomIdx() in left_atoms]
            # ��ʣ�ಿ��ת��ΪSMILES�ַ���
            left_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(left_atoms), bondsToUse=left_bonds)
            left_mol = Chem.MolFromSmiles(left_smiles, sanitize=False)
            left_smiles = Chem.MolToSmiles(left_mol)
            if '.' not in left_smiles:
                #��ɾ������ת��Ϊmol����
                groups_mol = get_groups_mol(mol, left_atoms)

                #print(f"Iteration {iteration}: {left_smiles}")
                mol = Chem.MolFromSmiles(left_smiles, sanitize=False)
                
                if mol is None:
                    raise ValueError(f"Failed to create molecule from SMILES: {left_smiles}")
                    
                #mol,_ = process_molecule(mol)
                for atom in mol.GetAtoms():
                    if atom.GetAtomMapNum() == 1:
                        atom_need = atom
                        #print(Chem.MolToSmiles(mol))
                        atom.SetAtomMapNum(0)
                        
                        info = get_atom_info(atom_need)
                        #print(info)
                        mol,test_smiles = process_molecule(mol)
                        #print(test_smiles)
                        '''
                        test_mol,_ = process_molecule(mol)
                        for i, atom in enumerate(test_mol.GetAtoms()):
                            print(bfs(test_mol, atom))
                            atom.SetAtomMapNum(i + 1)
                        
                        test_smiles = Chem.MolToSmiles(test_mol)
                        print(test_smiles)
                        '''
                        mol,_ = process_molecule(mol)
                        for atom_i in mol.GetAtoms():
                            #print(get_atom_info(atom_i))
                            if get_atom_info(atom_i) == info:
                                connected_atoms.append(atom_i.GetIdx())
                                break
                                
                #groups_mol,_ = process_molecule(groups_mol)
                for atom in groups_mol.GetAtoms():
                    if atom.GetAtomMapNum() == 2:
                        atom_need = atom
                        #print(Chem.MolToSmiles(groups_mol))
                        atom.SetAtomMapNum(0)
                        
                        groups_mol,test_smiles = process_molecule(groups_mol)
                        #print(test_smiles)
                        
                        info = get_atom_info(atom_need)
                        #print(info)
                        groups_mol,_ = process_molecule(groups_mol)
                        for atom_i in groups_mol.GetAtoms():
                            #print(get_atom_info(atom_i))
                            if get_atom_info(atom_i) == info:
                                connected_group_atoms.append(atom_i.GetIdx())
                                break

                successful_groups.append((groups_mol, connected_atoms, connected_group_atoms))
                atom_groups = get_atom_groups(mol)
                
                
                removed = True
                break
        if not removed:
            break
    return successful_groups, start_mol



def add_group(mol, group, connected_atoms, connected_group_atoms, bond_type):
    # �ϲ���������
    new_mol = Chem.CombineMols(mol, group)
    
    # ����һ��EditableMol����
    em = Chem.EditableMol(new_mol)
    
    # ȷ����������
    if bond_type == 1:
        bond_order = Chem.BondType.SINGLE
    elif bond_type == 2:
        bond_order = Chem.BondType.DOUBLE
    elif bond_type == 3:
        bond_order = Chem.BondType.TRIPLE
    else:
        raise ValueError(f"Invalid bond type: {bond_type}")
    
    # ������ԭ��������з���
    offset = mol.GetNumAtoms()
    for atom_idx, group_atom_idx in zip(connected_atoms, connected_group_atoms):
        em.AddBond(atom_idx, group_atom_idx+offset, order=bond_order)
    
    # ��ȡ�·���
    new_mol = em.GetMol()
    
    return new_mol

def replace_smiles(smiles):
    end_smiles = smiles
    end_smiles = end_smiles.replace("[CH3]", "C").replace("[CH2]", "C").replace("[CH]", "C").replace("[C]", "C")
    end_smiles = end_smiles.replace("[NH4]", "N").replace("[NH3]", "N").replace("[NH2]", "N").replace("[NH]", "N").replace("[N]", "N")
    end_smiles = end_smiles.replace("[SH5]", "S").replace("[SH4]", "S").replace("[SH3]", "S").replace("[SH2]", "S").replace("[SH]", "S").replace("[S]", "S")
    smiles = end_smiles.replace("[PH5]", "P").replace("[PH4]", "P").replace("[PH3]", "P").replace("[PH2]", "P").replace("[PH]", "P").replace("[P]", "P")
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    smiles = Chem.MolToSmiles(mol)
    return smiles

def is_valid_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        new_smiles = Chem.MolToSmiles(mol)
        return True
    except:
        return False





