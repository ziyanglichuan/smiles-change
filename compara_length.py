# -*- coding:gbk -*-
from smiles_change import *
import random

def read_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    data = {}
    for line in lines:
        smiles, code = line.strip().split()
        data[smiles] = code
    return data

vocab_data = read_file('./vocab/vocab_2.txt')


with open('./data/test_2.txt', 'r') as file:
    smiles_list = file.read().split('\n')


total_reduction = 0
num_pairs = 0

with open('./encode_data/out_code.txt', 'w') as outfile:
    smiles_list = random.sample(smiles_list, 10000)
    for smiles in smiles_list:
        if not is_valid_smiles(smiles) or smiles == "":
            print(f"Invalid SMILES: {smiles}")
            continue

        # 从SMILES字符串创建分子
        canonical_smiles = replace_smiles(smiles)
        mol = Chem.MolFromSmiles(canonical_smiles)

        atom_groups = get_atom_groups(mol)

        success = False
        for i in range(2):
            # 调用remove_groups函数
            removed_groups, start_groups = remove_groups(mol, atom_groups)

            # 逆反remove_groups函数的过程
            new_mol = start_groups

            for group, connected_atoms, connected_group_atoms in reversed(removed_groups):
                new_mol = add_group(new_mol, group, connected_atoms, connected_group_atoms, bond_type=1)
                new_smiles = Chem.MolToSmiles(new_mol)
                new_mol = Chem.MolFromSmiles(new_smiles, sanitize=False)

            end_smiles = Chem.MolToSmiles(new_mol)
            end_smiles = replace_smiles(end_smiles)

            if canonical_smiles == end_smiles:
                success = True
                break
            else:
                atom_groups.reverse()

        if not success:
            print(f"Failed case: original SMILES: {canonical_smiles}, transformed SMILES: {end_smiles}")
        else:
            len_result = 0
            result_str = ''
            start_smiles = Chem.MolToSmiles(start_groups)
            if start_smiles in vocab_data:
                code = vocab_data[start_smiles]
                result_str += code
                len_result += 1
            else:
                print(f'{start_smiles} not found in vocab_file')
            for group, connected_atoms, connected_group_atoms in reversed(removed_groups):
                group_smiles = Chem.MolToSmiles(group)
                if group_smiles in vocab_data:
                    code = vocab_data[group_smiles]
                    result_str += f'{connected_atoms[0]}/{connected_group_atoms[0]}'
                    result_str += code
                    len_result += 4
                else:
                    print(f'{group_smiles} not found in vocab_file')

            len_canonical = len(canonical_smiles)
            reduction = (len_canonical - len_result) / len_canonical * 100
            total_reduction += reduction
            num_pairs += 1

            outfile.write(f"{canonical_smiles} {result_str}\n")

average_reduction = total_reduction / num_pairs
print(f'Average reduction: {average_reduction:.2f}%')

