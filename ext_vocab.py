# -*- coding: gbk -*-
from smiles_change import *
import argparse

def extVocab_smiles(smiles_list):
    vocab_set = set()
    for i, smiles in enumerate(smiles_list):
        if not is_valid_smiles(smiles) or smiles == "":
            print(f"Invalid SMILES: {smiles}")
            continue

        canonical_smiles = replace_smiles(smiles)
        mol = Chem.MolFromSmiles(canonical_smiles)

        try:
            atom_groups = get_atom_groups(mol)

            # 调用remove_groups函数
            removed_groups, start_groups = remove_groups(mol, atom_groups)

        except Exception as e:
            print(f"Error processing SMILES at line {i+1}: {smiles}")
            print(f"Error: {e}")
            continue
            
        start_smiles = Chem.MolToSmiles(start_groups)
        vocab_set.add(start_smiles)
        
        # 打印结果
        for group, connected_atoms, connected_group_atoms in removed_groups:
            group_smiles = Chem.MolToSmiles(group)
            vocab_set.add(group_smiles)

    vocab_list = sorted(list(vocab_set))
    return vocab_list

def encode_vocab(vocab_list):
    vocab_dict = {}
    chars = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
    for i, code in enumerate(itertools.product(chars, repeat=3)):
        if i >= len(vocab_list):
            break
        vocab_dict[vocab_list[i]] = ''.join(code)
    return vocab_dict


def ext_file(input_filename, output_filename):
    vocab_set = set()
    with open(input_filename, 'r') as file:
        smiles_list = file.read().split('\n')

    vocab_list = extVocab_smiles(smiles_list)
    vocab_dict = encode_vocab(vocab_list)

    with open(output_filename, 'w') as outfile:
        for smiles, code in vocab_dict.items():
            outfile.write(f'{smiles} {code}\n')
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_filename', help='Input file name')
    parser.add_argument('output_filename', help='Output file name')
    args = parser.parse_args()

    ext_file(args.input_filename, args.output_filename)
