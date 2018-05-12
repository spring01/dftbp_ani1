"""
Usage: python check_agreement.py gdict_dftb.p ./mio-0-1/ ./dftbplus-18.1.x86_64-linux/bin/dftb+
"""

import os
import sys
import subprocess
import re
import collections
import cPickle as pickle


def calc_dftbplus(geom, bin_dftbplus):
    """Will run the calculation and write log files to the current path"""
    # write dftb_in.hsd input file
    write_input(geom)

    # run dftbplus
    proc = subprocess.Popen([bin_dftbplus], stdout=subprocess.PIPE)
    output, err = proc.communicate()

    # read and parse output
    with open('detailed.out', 'r') as out:
        for line in out:
            if 'Total Electronic energy' in line:
                elec = float(re.findall('[+-]?\d+\.\d+', line)[0])
            if 'Repulsive energy' in line:
                rep = float(re.findall('[+-]?\d+\.\d+', line)[0])
    return elec, rep

def write_input(geom):
    """Write dftb_in.hsd to the current path"""
    z2atom = {1: 'H', 6: 'C', 7: 'N', 8: 'O'}
    atom2ang = {'H': 'H = "s"', 'C': 'C = "p"', 'N': 'N = "p"', 'O': 'O = "p"'}
    atom_list = [z2atom[charge] for charge in geom.z]
    atom_set = set(atom_list)
    atom_ulist = sorted(list(atom_set))
    skfile_list = []
    ang_list = []
    for a1 in atom_ulist:
        ang_list.append(atom2ang[a1])
        for a2 in atom_ulist:
            skname = '{}-{}'.format(a1, a2)
            skfile_path = os.path.join(param_dir, '{}.skf'.format(skname))
            skfile_list.append('{} = "{}"'.format(skname, skfile_path))

    atom_counter = sorted(collections.Counter(atom_list).items())
    atom2idx = {atom: idx + 1 for idx, (atom, _) in enumerate(atom_counter)}
    geom_list = []
    for iatom, (charge, xyz) in enumerate(zip(geom.z, geom.rcart.T)):
        atom = z2atom[charge]
        coords = '{:18.11E} {:18.11E} {:18.11E}'.format(*(xyz))
        geom_list.append('{} {} {}'.format(iatom + 1, atom2idx[atom], coords))

    # Geometry block
    block_geom = 'Geometry = GenFormat {\n'
    block_geom += '{} C\n '.format(len(atom_list))
    for atom in atom_ulist:
        block_geom += ' {}'.format(atom)
    block_geom += '\n\n'
    for atom_coords in geom_list:
        block_geom += '  {}\n'.format(atom_coords)
    block_geom += '}'

    # Hamiltonian block
    block_ham = 'Hamiltonian = DFTB {\n  SCC = Yes\n  SlaterKosterFiles {\n'
    for skfile in skfile_list:
        block_ham += '    {}\n'.format(skfile)
    block_ham += '  }\n  MaxAngularMomentum {\n'
    for ang in ang_list:
        block_ham += '    {}\n'.format(ang)
    block_ham += '  }  OrbitalResolvedSCC = Yes\n'
    block_ham += '\n}\n\nParserOptions {\n  ParserVersion = 5\n}\n'

    # write input
    dftb_in = '{}\n\n{}'.format(block_geom, block_ham)
    with open('dftb_in.hsd', 'w') as hsd:
        hsd.write(dftb_in)




if __name__ == '__main__':
    gdict_filename = sys.argv[1]
    param_dir = sys.argv[2]
    bin_dftbplus = sys.argv[3]

    with open(gdict_filename,'rb') as f:
        gdict_dftb = pickle.load(f)

    print '%s, %s, %s, %s, %s, %s, %s' % ('dataset', 'ibatch', 'imol', 'elec', 'elec_dftbp', 'rep', 'rep_dftbp')
    for dataset in 'train', 'test':
        for ibatch in range(200):
            for imol in range(66):
                mol = gdict_dftb[dataset][ibatch][imol]
                elec = mol['Edftb_elec']
                rep = mol['Edftb_rep']
                geom = mol['geom']
                ref_elec, ref_rep = calc_dftbplus(geom, bin_dftbplus)
                print '%s, %3d, %2d, %14.8f, %14.8f, %14.8f, %14.8f' % (dataset, ibatch, imol, elec, ref_elec, rep, ref_rep)


