# coding: utf-8

import MDAnalysis
import sys

sys.setrecursionlimit(100000)



def parser():
    list_name  = sys.argv[1]
    pdb_name = sys.argv[2]

    if len(sys.argv) == 4:
        out_name = sys.argv[3]
    else :
        out_name = list_name

    return list_name, pdb_name, out_name

def main(list_name, pdb_name, out_name):

    atoms_list = open(list_name, 'r').readlines()

    atoms_num = []

    for i in range(len(atoms_list)):
        atoms_ar = atoms_list[i].split()

        for j in range(len(atoms_ar)):
            try :
                atoms_num.append(int(atoms_ar[j]))
            except ValueError:
                pass

    print('The zone has %i atoms (if it is a QM zone, link atoms are not counted)' % len(atoms_num))

    u = MDAnalysis.Universe(pdb_name)

    atoms_sel_list = 'bynum ' + ' or bynum '.join([ str(at) for at in atoms_num ])


    #for i in range(len(atoms_num)):
    #    if i == len(atoms_num) - 1:
    #        atoms_sel_list = atoms_sel_list + 'bynum ' + str(atoms_num[i])
    #    else :
    #        atoms_sel_list = atoms_sel_list + 'bynum ' + str(atoms_num[i]) + ' or '

    atoms_sel = u.select_atoms(atoms_sel_list)

    atoms_sel.write(out_name + '.pdb')


if __name__ == '__main__':
    list_name, pdb_name, out_name = parser()
    main(list_name, pdb_name, out_name)

# usage:
# conda activate md_traj (venv with MDAnalysis)
# python atoms_list_to_pdb.py 'qm_list' open_hf_11677_solv.cropped.pdb qm_region
