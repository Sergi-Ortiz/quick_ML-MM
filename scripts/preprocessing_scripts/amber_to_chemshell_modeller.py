#!/home/mcanyelles/miniconda3/envs/pyenv/bin/python
# -*- coding: utf-8 -*-

# TODO check Fe and O

# Import packages
from parmed import load_file
from MDAnalysis import Universe
from MDAnalysis.exceptions import SelectionError
import sys, os
from argparse import ArgumentParser

def parser():
    parser = ArgumentParser(description="Script to adapt AMBER topology and coordinates to be used in ChemShell.")

    parser.add_argument('-p', '--parameters',
                        help='AMBER prmtop file containing topology and parameters of the system.',
                        type=str,
                        required=True
                        )

    parser.add_argument('-c', '--coordinates',
                        help='AMBER inpcrd file containing the coordinates of the system',
                        required=True,
                        type=str
                        )


    parser.add_argument('-cs', '--crop_system',
                        help="Trigger for activating the cropping of the system's solvent in a radius around a residue. If used, '-csc' has to be specified and '-csr' is optional with a default value of 17 A",
                        action='store_true',
                        default=False,
                        required=False
                        )
    parser.add_argument('-csc', '--crop_system_center',
                        help='Number of the residue that will be used for creating the drop of solvent',
                        required=False,
                        default=0
                        )
    parser.add_argument('-csr', '--crop_system_radius',
                        help="Radius (in Å) around the central residue specified with the '-csc' flap. Default value is 17 Å.",
                        required=False,
                        default=17,
                        type=float
                        )

    parser.add_argument('-al', '--active_atoms_list',
                        help="Trigger for creating a tcl list containing the list of atoms around an specific atom. This is required for ChemShell QM/MM simulations.",
                        required=False,
                        action='store_true',
                        default=False
                        )
    parser.add_argument('-alc', '--active_atoms_list_center',
                        help='Number of the residue that will be used as centre for creating the active atoms list',
                        required=False,
                        default=0,
                        )
    parser.add_argument('-alr', '--active_atoms_list_radius',
                        help="Radius (in Å) around the central atom specified with the '-alc' flap. Default value is 15 Å.",
                        required=False,
                        default=15,
                        type=float
                        )

    parser.add_argument('-ra', '--rename_atoms',
                        help="Trigger for changing conflicting atom names/types used in AMBER to atoms names understandable by ChemShell. Atom names from MCPB.py are also renamed if the original names are given with the '-ram' flag.",
                        required=False,
                        action='store_true',
                        default=False
                        )
    parser.add_argument('-ram', '--rename_atoms_MCPB',
                        help="List of atom names for changing MCPB names (like Y1, Y2, M1, etc.)",
                        nargs='+',
                        type=str,
                        default=None,
                        required=False
                        )
    parser.add_argument('-rac', '--rename_atoms_MCPB_checker',
                        help="Trigger for checking the existance of any atom name unknown by ChemShell. It prints the list of unknown atoms and exists the program.",
                        required=False,
                        default=False,
                        action='store_true'
                        )

    parser.add_argument('-o', '--output',
                        help='Name of the output files.',
                        type=str,
                        default=None
                        )


    args = parser.parse_args()

    if args.crop_system == True and args.crop_system_center == 0:
        sys.exit("If '-cs' is activated, the center residue for cropping the solvent has to be specified with the '-csc' flag.")

    if args.active_atoms_list == True and args.active_atoms_list_center == 0:
        sys.exit("If '-al' is activated, the center atom for creating the atoms list has to be specified with the '-alc' flag.")

    if args.rename_atoms and args.rename_atoms_MCPB == None:
        print('MCPB atom types will be automatically adapted to ChemShell atom types.')

    if args.output == None:
        args.output = args.parameters[:-7]

    return args



def crop_topology(args):

    # load paramters
    u_top = Universe(args.parameters)

    #
    try :
        sel = str(u_top.select_atoms('resid %s' % args.crop_system_center).residues)

    except SelectionError:
        sys.exit("The selected residue for cropping the topology does not exist.")


    del u_top, sel

    topology = load_file(args.parameters, args.coordinates)

    topology.box = None
    topology.strip(':WAT&!:%s<:%s' % (args.crop_system_center, args.crop_system_radius))
    topology.strip(':Na+&!:%s<:%s' % (args.crop_system_center, args.crop_system_radius))
    topology.strip(':Cl-&!:%s<:%s' % (args.crop_system_center, args.crop_system_radius))
    topology.strip(':K+&!:%s<:%s' % (args.crop_system_center, args.crop_system_radius))

    topology.write_parm(args.output + '_cropped.prmtop')
    topology.save(args.output + '_cropped.inpcrd')
    #for res in topology.residues:
    #   res.ter = False
    topology.write_pdb(args.output + '_cropped.pdb')

    pdb = open(args.output + '_cropped.pdb', 'r').readlines()
    pdb_o = open(args.output + '_cropped.pdb.tmp', 'w')

    atom = int(pdb[0][6:11])
    for line in pdb:
        if 'TER' in line:
            pass
        else :
            new_line = line[0:6] + '{:5}' + line[11:]
            pdb_o.write(new_line.format(atom))
            atom += 1

    pdb_o.close()
    os.remove(args.output + '_cropped.pdb')
    os.rename(args.output + '_cropped.pdb.tmp', args.output + '_cropped.pdb')

    print('Cropped topology and coordinates have been saved as \'%s\', \'%s\' and \'%s\'' % (args.output + '_cropped.prmtop', args.output + '_cropped.inpcrd', args.output + '_cropped.pdb'))

    args.parameters = args.output + '_cropped.prmtop'
    args.coordinates = args.output + '_cropped.inpcrd'
    return args


def create_active_atoms_list(args):

    u_set_act = Universe(args.parameters, args.coordinates)

    try :
        sel_carbon = str(u_set_act.select_atoms("bynum %s" % args.active_atoms_list_center))
        locA = sel_carbon.find('[<') + 2
        locB = sel_carbon.find(' and segid')

    except SelectionError:
        sys.exit('The atom center for creating the active atoms list does not exists.')

    selection = u_set_act.select_atoms(str('byres around %s bynum %s' % (args.active_atoms_list_radius, args.active_atoms_list_center)))

    txt = open(args.output + '.act_list', 'w')
    txt.write('set act { ')
    for i in range(0, len(selection)):
        locA = str(selection[i]).find('<Atom ') + 6
        locB = str(selection[i]).find(': ')
        txt.write(str(selection[i])[locA:locB] + " ")
    txt.write("} ")
    txt.close()

    print('%s atoms have been set as active for the QM/MM calculation using ChemShell.' % (len(selection)))


def topology_checker(args):

    top    = open(args.parameters, 'r').readlines()
    top_out = open(str(args.parameters)[:-7] + '_mod.prmtop', 'w')

    initial = 0
    heterotypes_in = []
    metals_in      = []
    for i in range(0, len(top)):
        if '%FLAG AMBER_ATOM_TYPE' in top[i]:
            initial = i +1
        if '%FLAG TREE_CHAIN_CLASSIFICATION' in top[i]:
            final   = i
        if ('  Y' in top[i] or top[i].find('Y') == 0) and '%' not in top[i]:
            index = 0
            while index < len(top[i]):
                index = top[i].find('Y', index)
                if index == -1:
                    break
                elif index != -1:
                    if top[i].find('Y') == 0:
                        loc = 0
                        heterotypes_in.append(str(top[i])[loc:loc+3])
                    else :
                        loc = top[i].find(' Y', index-1) + 1
                        #print(loc)
                        heterotypes_in.append(str(top[i])[loc:loc+3])

                    index += 1

        if ' M' in top[i] and '%' not in top[i]:
            loc = top[i].find(' M') + 1
            try :
                int(str(top[i])[loc+1:loc+3])
                metals_in.append(str(top[i])[loc:loc+3])
            except ValueError:
                pass


def topology_adapter(args):

    new_types = []

    u = Universe(args.parameters)

    top    = open(args.parameters, 'r').readlines()
    top_out = open(str(args.parameters)[:-7] + '_mod.prmtop', 'w')

    initial = 0
    heterotypes_in = []
    metals_in      = []
    for i in range(0, len(top)):
        if '%FLAG AMBER_ATOM_TYPE' in top[i]:
            initial = i +1
        if '%FLAG TREE_CHAIN_CLASSIFICATION' in top[i]:
            final   = i
        if ('  Y' in top[i] or top[i].find('Y') == 0) and '%' not in top[i]:
            index = 0
            while index < len(top[i]):
                index = top[i].find('Y', index)
                if index == -1:
                    break
                elif index != -1:
                    if top[i].find('Y') == 0:
                        loc = 0
                        heterotypes_in.append(str(top[i])[loc:loc+3])
                    else :
                        loc = top[i].find(' Y', index-1) + 1
                        #print(loc)
                        heterotypes_in.append(str(top[i])[loc:loc+3])

                    index += 1

        if ' M' in top[i] and '%' not in top[i]:
            loc = top[i].find(' M') + 1
            try :
                int(str(top[i])[loc+1:loc+3])
                metals_in.append(str(top[i])[loc:loc+3])
            except ValueError:
                pass
        
    #print(heterotypes_in)
    #print(metals_in)

    if args.rename_atoms_MCPB == None and (len(heterotypes_in) > 0 or len(metals_in) > 0):

        translator = {
            'NE'  : 'N ',
            'NZ'  : 'N ',
            'NE1' : 'N  ',
            'NE2' : 'N  ',
            'ND1' : 'N  ',
            'ND2' : 'N  ',
            'O'   : 'O',
            'OW'  : 'O ',
            'OD1' : 'O  ',
            'OD2' : 'O  ',
            'OE1' : 'O  ',
            'OE2' : 'O  ',
            'OXT' : 'O  ',
            'SG'  : 'S ',
            'SD'  : 'S ',
        }

        args.rename_atoms_MCPB = []
        for t in heterotypes_in:
     #       print(u.select_atoms(f"type {t}"))
            args.rename_atoms_MCPB.append(translator[u.select_atoms(f"type {t}").names[0]])

        for t in metals_in:
            args.rename_atoms_MCPB.append(u.select_atoms(f"type {t}").names[0])


    for l in range(0, initial):
        top_out.write(top[l])

    for l in range(initial, final):
        l_ = top[l].replace('hc', 'H ')
        l_ = l_.replace('ha', 'H ')
        l_ = l_.replace('h1', 'H ')
        l_ = l_.replace('ho', 'H ')

        l_ = l_.replace('2C', 'C2')
        l_ = l_.replace('3C', 'C3')

        l_ = l_.replace('CO', 'C ')
        l_ = l_.replace('CX', 'C ')
        l_ = l_.replace('c2', 'C ')
        l_ = l_.replace('c3', 'C ')
        l_ = l_.replace('cx', 'C ')
        l_ = l_.replace('ce', 'C ')
        l_ = l_.replace('cf', 'C ')

        l_ = l_.replace('op', 'O ')
        l_ = l_.replace('os', 'O ')
        l_ = l_.replace('oh', 'O ')

        l_ = l_.replace('c ', 'C ')
        l_ = l_.replace('o ', 'O ')

        l_ = l_.replace('Na+', 'NA+')
        l_ = l_.replace('Cl-', 'CL-')
        l_ = l_.replace('K+', 'K+')
        l_ = l_.replace('Ca+', 'Ca+')

        for j in range(len(args.rename_atoms_MCPB)):
            whitespaces = (len((heterotypes_in + metals_in)[j]) - len(args.rename_atoms_MCPB[j]))*' '
            l_ = l_.replace((heterotypes_in + metals_in)[j], args.rename_atoms_MCPB[j] + whitespaces)

        top_out.write(l_)

    for l in range(final, len(top)):
        top_out.write(top[l])

    top_out.close()

    if args.rename_atoms_MCPB != None:
        print(f"{len(args.rename_atoms_MCPB)} atom types from MCPB have been changed. ")



def main():
    args = parser()

    if args.crop_system:
        args = crop_topology(args)
        print('New parameters and coordinates files of the cropped system will be used in future actions.')

    if args.active_atoms_list:
        create_active_atoms_list(args)

    if args.rename_atoms:
        topology_adapter(args)
        print('prmtop fixed!')




if __name__ == '__main__':
    main()
