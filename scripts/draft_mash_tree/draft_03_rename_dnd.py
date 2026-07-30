'''
Rename accession numbers in the mashtree output to sample names.
'''

import os

def rename_tree():
    mapping_file = os.environ.get("MAPPING_FILE", "./datasets/name_mapping.csv")
    tree_file    = os.environ.get("TREE_FILE", "./mashtree.dnd")
    output_file  = os.environ.get("OUTPUT_FILE", "./renamed_tree.dnd")

    replacements = []
    with open(mapping_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and ',' in line:
                srr, name = line.split(',')
                replacements.append((srr.strip(), name.strip()))

    with open(tree_file, 'r') as f:
        tree_content = f.read()

    for srr, name in replacements:
        tree_content = tree_content.replace(srr, name)
        tree_content = tree_content.replace(".unitigs", "")

    with open(output_file, 'w') as f:
        f.write(tree_content)

    print(f"Done. Wrote {output_file}")


if __name__ == "__main__":
    rename_tree()