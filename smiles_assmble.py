# -*- coding:gbk -*-
from smiles_change import *


# ��SMILES�ַ�����������
smiles = "CCCC1CCC(C[NH3+])C([NH+]2CC(C)CC(C)C2)C1"

canonical_smiles = replace_smiles(smiles)
mol = Chem.MolFromSmiles(canonical_smiles)
print("smiles:",canonical_smiles)

atom_groups = get_atom_groups(mol)
#print(atom_groups)
#���ֻ���ת��ΪSMILES�ַ���
for groups in atom_groups:
    smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(groups))
    print("atom_groups:",smiles)


#atom_groups.reverse()

# ����remove_groups����
removed_groups, start_groups = remove_groups(mol, atom_groups)

# ��ӡ���
for group, connected_atoms, connected_group_atoms in removed_groups:
    group_smiles = Chem.MolToSmiles(group)
    print(f"Removed group: {group_smiles}, connected to: {connected_atoms},group_connect:{connected_group_atoms}")

# �淴remove_groups�����Ĺ���
new_mol = start_groups
new_smiles = Chem.MolToSmiles(new_mol)
print(new_smiles)

for group, connected_atoms, connected_group_atoms in reversed(removed_groups):
   
    # ���connected_atoms��Ϊ�գ���ԭ���鰴��connected_atoms��λ�����ӵ���ʼ�����������·���
    new_mol = add_group(new_mol, group, connected_atoms, connected_group_atoms,bond_type=1)
    
    # ��ȡnew_mol���ӵ�SMILES�ַ���
    new_smiles = Chem.MolToSmiles(new_mol)
    new_mol = Chem.MolFromSmiles(new_smiles, sanitize=False)

    new_smiles = Chem.MolToSmiles(new_mol)
    print(new_smiles)
    
end_smiles = replace_smiles(new_smiles)
print(end_smiles)

with open('./data/test.txt', 'r') as file:
    smiles_list = file.read().split('\n')


with open('out.txt', 'w') as outfile:
    success_count = 0
    counter = 0
    for smiles in smiles_list:
    
        if not is_valid_smiles(smiles) or smiles == "":
            print(f"Invalid SMILES: {smiles}")
            continue
            
        counter += 1
        # ��SMILES�ַ�����������
        canonical_smiles = replace_smiles(smiles)
        mol = Chem.MolFromSmiles(canonical_smiles)
        
        atom_groups = get_atom_groups(mol)
        
        success = False
        for i in range(2):
            # ����remove_groups����
            removed_groups, start_groups = remove_groups(mol, atom_groups)
            
            # �淴remove_groups�����Ĺ���
            new_mol = start_groups
            
            for group, connected_atoms, connected_group_atoms in reversed(removed_groups):
                new_mol = add_group(new_mol, group, connected_atoms, connected_group_atoms, bond_type=1)
                new_smiles = Chem.MolToSmiles(new_mol)
                new_mol = Chem.MolFromSmiles(new_smiles, sanitize=False)
            
            end_smiles = Chem.MolToSmiles(new_mol)
            end_smiles = replace_smiles(end_smiles)
            
            if canonical_smiles == end_smiles:
                success_count += 1
                success = True
                break
            else:
                atom_groups.reverse()
        
        if not success:
            print(f"Failed case: original SMILES: {canonical_smiles}, transformed SMILES: {end_smiles}")
    
        if counter % 5000 == 0:
            success_rate = success_count / counter
            outfile.write(f"Success rate after processing {counter} SMILES: {success_rate}\n")
            outfile.flush()
    
    success_rate = success_count / len(smiles_list)
    print(f"Final success rate: {success_rate}")