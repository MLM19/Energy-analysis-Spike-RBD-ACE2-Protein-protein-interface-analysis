import argparse
import os
import sys
import math
import matplotlib.pyplot as plt

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.PDBIO import PDBIO, Select
from Bio.PDB.NeighborSearch import NeighborSearch

# set the pdb_path and load the structure
pdb_path = "/Users/msteph.padilla/Desktop/BDBI/second_year/First/Biophysics/Project/7wbq_fixed.pdb"

class ResiduesDataLib():
    def __init__(self, fname):
        self.residue_data = {}
        try:
            fh = open(fname, "r")
        except OSError:
            print("#ERROR while loading library file (", fname, ")")
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
            r = Residue(data)
            self.residue_data[r.id] = r
        self.nres = len(self.residue_data)

    def get_params(self, resid, atid):
        atom_id = resid + ':' + atid
        if atom_id in self.residue_data:
            return self.residue_data[atom_id]
        else:
            print("WARNING: atom not found in library (", atom_id, ')')
            return None

class Residue():
    def __init__(self,data):
        self.id     = data[0]+':'+data[1]
        self.at_type = data[2]
        self.charge  = float(data[3])
        
class AtomType():
    def __init__(self, data):
        self.id   = data[0]
        self.eps  = float(data[1])
        self.sig  = float(data[2])
        self.mass = float(data[3])
        self.fsrf = float(data[4])
        self.rvdw = self.sig * 0.5612
        
class VdwParamset(): #extracted from GELPI's github
    #parameters for the VdW
    def __init__ (self, file_name):
        self.at_types = {}
        try:
            fh = open(file_name, "r")
        except OSError:
            print ("#ERROR while loading parameter file (", file_name, ")")
            sys.exit(2)
        for line in fh:
            if line[0] == '#':
                continue
            data = line.split()
            self.at_types[data[0]] = AtomType(data)
        self.ntypes = len(self.at_types)
        fh.close()

NACCESS_BINARY = '/Users/msteph.padilla/Desktop/BDBI/second_year/First/Biophysics/Project/soft/NACCESS/naccess'

parse_cmd = argparse.ArgumentParser(
    prog='binding',
    description='binding energy calculation'
)

parse_cmd.add_argument(
    '--rlib',
    action='store',
    dest='reslib_file',
    default='data/aaLib.lib',
    help='Residue Library'
)
parse_cmd.add_argument(
    '--vdw',
    action='store',
    dest='vdwprm_file',
    default='data/vdwprm',
    help='Vdw parameters'
)

parse_cmd.add_argument(
    '--dist',
    action='store',
    dest='cutoff_dist',
    default=4.3,
    type=float,
    help='Cutoff distance for determining the interface (0: use all residues):'
)
parse_cmd.add_argument(
    'pdb_path',
    nargs='?',
    const="/Users/msteph.padilla/Desktop/BDBI/second_year/First/Biophysics/Project/6m0j_fixed.pdb",
    help='Input PDB'
)
# #This line is not needed, because the path is already defined.
    #parse_cmd.add_argument('pdb_path', help='Input PDB', type=open) 

args = parse_cmd.parse_args()


# loading residue library from data/aaLib.lib
residue_library = ResiduesDataLib('/Users/msteph.padilla/Desktop/BDBI/second_year/First/Biophysics/Project/aaLib.lib')

# loading VdW parameters
ff_params = VdwParamset('/Users/msteph.padilla/Desktop/BDBI/second_year/First/Biophysics/Project/vdwprm')

# Setting the Bio.PDB.Parser object
parser = PDBParser(PERMISSIVE=1)

# Loading structure
st = parser.get_structure('STR', pdb_path)

ala_atoms = {'N', 'H', 'CA', 'HA', 'C', 'O', 'CB', 'HB', 'HB1', 'HB2', 'HB3', 'HA1', 'HA2', 'HA3'}

def residue_id(res):
    '''Returns readable residue id'''
    return '{} {}{}'.format(res.get_resname(), res.get_parent().id, res.id[1])

def atom_id(at):
    '''Returns readable atom id'''
    return '{}.{}'.format(residue_id(at.get_parent()), at.id)

#Here all the functions needed that can be found in the energy.py
def add_atom_parameters(st, res_lib, ff_params):
    ''' Adds parameters from libraries to atom .xtra field
        For not recognized atoms, issues a warning and put default parameters
    '''
    for at in st.get_atoms():
        resname = at.get_parent().get_resname()
        params = res_lib.get_params(resname, at.id)
        if not params:
            at.xtra['atom_type'] = at.element
            at.xtra['charge'] = 0
        else:
            at.xtra['atom_type'] = params.at_type
            at.xtra['charge'] = params.charge
        at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]

def get_interface(st, dist):
    ''' Detects interface residues within a distance(dist)
        Assumes two chains, i.e. a unique interface set per chain.
    '''
    select_ats = []
    for at in st.get_atoms():
        # Skip Hydrogens to reduce time
        if at.element != 'H':
            select_ats.append(at)
    nbsearch = NeighborSearch(select_ats)
    interface = {}
    # Sets are more efficient than lists. Use sets when order is not relevant
    for ch in st[0]:
        interface[ch.id] = set()

    for at1, at2 in nbsearch.search_all(dist):
        #Only different chains
        res1 = at1.get_parent()
        ch1 = res1.get_parent()
        res2 = at2.get_parent()
        ch2 = res2.get_parent()
        if ch1 != ch2:
            interface[ch1.id].add(res1)
            interface[ch2.id].add(res2)
    return interface

def calc_solvation(st,res):
    '''Solvation energy based on ASA'''
    solv = 0.
    solv_ala = 0.
    for at in res.get_atoms():
        if 'EXP_NACCESS' not in at.xtra:
            continue
        s = float(at.xtra['EXP_NACCESS'])* at.xtra['vdw'].fsrf
        solv += s
        if at.id in ala_atoms:
            solv_ala += s
    return solv, solv_ala

def calc_int_energies(st, res):
    '''Returns interaction energies (residue against other chains)
        for all atoms and for Ala atoms
    '''
    elec = 0.
    elec_ala = 0.
    vdw = 0.
    vdw_ala = 0.

    for at1 in res.get_atoms():
        for at2 in st.get_atoms():
        # skip same chain atom pairs
            if at2.get_parent().get_parent() != res.get_parent(): ##2
                r = at1 - at2
                e = elec_int(at1, at2, r)
                elec += e
                if at1.id in ala_atoms: #GLY are included implicitly
                     elec_ala += e
                e = vdw_int(at1, at2, r)
                vdw += e
                if at1.id in ala_atoms: #GLY are included implicitly
                     vdw_ala += e
    return elec, elec_ala, vdw, vdw_ala

def MH_diel(r):
    '''Mehler-Solmajer dielectric'''
    return 86.9525 / (1 - 7.7839 * math.exp(-0.3153 * r)) - 8.5525

def elec_int(at1, at2, r):
    '''Electrostatic interaction energy between two atoms at r distance'''
    return 332.16 * at1.xtra['charge'] * at2.xtra['charge'] / MH_diel(r) / r

def vdw_int(at1, at2, r):
    '''Vdw interaction energy between two atoms'''
    eps12 = math.sqrt(at1.xtra['vdw'].eps * at2.xtra['vdw'].eps)
    sig12_2 = at1.xtra['vdw'].sig * at2.xtra['vdw'].sig
    return 4 * eps12 * (sig12_2**6/r**12 - sig12_2**3/r**6)

# assign data types, and charges from libraries
# We will use the xtra attribute in Bio.PDB.Atom to hold the new data
# Possible errors on N-term and C-Term atoms
# Possible errors on HIS alternative forms
add_atom_parameters(st, residue_library, ff_params )

# Calculating surfaces
# The specific PATH to naccess script (in soft) is needed
# ASA goes to .xtra field directly

srf = NACCESS_atomic(st[0], naccess_binary=NACCESS_BINARY)

# Prepare surfaces for the separate chains
# Alternatively the twp PDB files can be prepared outside and parsed here
io = PDBIO()
st_chains = {}

# Using BioIO trick (see tutorial) to select chains
class SelectChain(Select):
    def __init__(self, chid):
        self.id = chid

    def accept_chain(self, chain):
        if chain.id == self.id:
            return 1
        else:
            return 0

for ch in st[0]:
    io.set_structure(st)
    io.save('tmp.pdb', SelectChain(ch.id))
    st_chains[ch.id] = parser.get_structure('stA', 'tmp.pdb')
    add_atom_parameters(st_chains[ch.id], residue_library, ff_params)
    srfA = NACCESS_atomic(st_chains[ch.id][0], naccess_binary=NACCESS_BINARY)
os.remove('tmp.pdb')

## Initiatlize Energy aggregates
elec = {}
elec_ala = {}

vdw = {}
vdw_ala = {}

solvAB = {}
solvAB_ala = {}

solvA = {}
solvA_ala = {}

totalIntElec = 0.
totalIntVdw = 0.
totalSolv = 0.
totalSolvMon = {}
## We get the chsin ids,not always they are A and B
chids = []
for ch in st[0]:
    chids.append(ch.id)
    totalSolvMon[ch.id] = 0

total = 0.

for ch in st[0]:
    for res in ch.get_residues():
        elec[res], elec_ala[res], vdw[res], vdw_ala[res] = calc_int_energies(st[0], res)
        solvAB[res], solvAB_ala[res] = calc_solvation(st[0], res)
        solvA[res], solvA_ala[res] = calc_solvation(
            st_chains[ch.id],
            st_chains[ch.id][0][ch.id][res.id[1]]
        )
        totalIntElec += elec[res]
        totalIntVdw += vdw[res]
        totalSolv += solvAB[res]
        totalSolvMon[ch.id] += solvA[res]
        total += elec[res] + vdw[res] + solvAB[res] - solvA[res]
print("Interaction energy for all residues")
print('{:20}: {:11.4f}'.format('Total Elec Int.', totalIntElec))
print('{:20}: {:11.4f}'.format('Total Vdw Int.', totalIntVdw))
print('{:20}: {:11.4f}'.format('Total Solv AB', totalSolv))
print('{:19}{}: {:11.4f}'.format('Total Solv ', chids[0], totalSolvMon[chids[0]]))
print('{:19}{}: {:11.4f}'.format('Total Solv ', chids[1], totalSolvMon[chids[1]]))
print('{:20}: {:11.4f}'.format('DGintAB-A-B', total))
print("\#R {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:11.4f}".format(totalIntElec, totalIntVdw, totalSolv, totalSolvMon[chids[0]], totalSolvMon[chids[1]], total))

#Calculation of interaction energies based in interface
if args.cutoff_dist > 0.:
    interface = get_interface(st, args.cutoff_dist)

## Initiatlize Energy aggregates
Intelec = {}
Intelec_ala = {}

Intvdw = {}
Intvdw_ala = {}

IntsolvAB = {}
IntsolvAB_ala = {}

IntsolvA = {}
IntsolvA_ala = {}

totalIntElec = 0.
totalIntVdw = 0.
totalSolv = 0.
totalSolvMon = {}

## We get the chsin ids,not always they are A and B
chids = []
for ch in st[0]:
    chids.append(ch.id)
    totalSolvMon[ch.id] = 0

total = 0.

for ch in st[0]:
    for res in ch.get_residues():
        if args.cutoff_dist > 0 and res not in interface[ch.id]:
            continue
        Intelec[res], Intelec_ala[res], Intvdw[res], Intvdw_ala[res] = calc_int_energies(st[0], res)
        IntsolvAB[res], IntsolvAB_ala[res] = calc_solvation(st[0], res)
        IntsolvA[res], IntsolvA_ala[res] = calc_solvation(
            st_chains[ch.id],
            st_chains[ch.id][0][ch.id][res.id[1]]
        )
        totalIntElec += Intelec[res]
        totalIntVdw += Intvdw[res]
        totalSolv += IntsolvAB[res]
        totalSolvMon[ch.id] += IntsolvA[res]
        total += Intelec[res] + Intvdw[res] + IntsolvAB[res] - IntsolvA[res]
print("Interaction energy based in interface residues only")
print('{:20}: {:11.4f}'.format('Total Elec Int.', totalIntElec))
print('{:20}: {:11.4f}'.format('Total Vdw Int.', totalIntVdw))
print('{:20}: {:11.4f}'.format('Total Solv AB', totalSolv))
print('{:19}{}: {:11.4f}'.format('Total Solv ', chids[0], totalSolvMon[chids[0]]))
print('{:19}{}: {:11.4f}'.format('Total Solv ', chids[1], totalSolvMon[chids[1]]))
print('{:20}: {:11.4f}'.format('DGintAB-A-B', total))
print("\#R {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:11.4f} {:11.4f}".format(totalIntElec, totalIntVdw, totalSolv, totalSolvMon[chids[0]], totalSolvMon[chids[1]], total))

# Store the calculated values
residue_ids = []
bar_values = []

print("Ala Scanning: DDGs for X->Ala mutations on interface residues")
for ch in st[0]:
    for resi in ch.get_residues():
        if args.cutoff_dist > 0 and resi not in interface[ch.id]:
            continue
        residue_id_str = residue_id(resi)
        print(
            '{:11} {:11.4f}{:11.4f}{:11.4f}{:11.4f}{:11.4f}'.format(
                residue_id(resi),
                - Intelec[resi] + elec_ala[resi],
                - Intvdw[resi] + Intvdw_ala[resi],
                - IntsolvAB[resi] + IntsolvAB_ala[resi],
                - IntsolvA[resi] + IntsolvA_ala[resi],
                - Intelec[resi] + Intelec_ala[resi] - Intvdw[resi] + Intvdw_ala[resi] - IntsolvAB[resi] +\
                    IntsolvAB_ala[resi] - IntsolvA[resi] + IntsolvA_ala[resi]
            )
        )
        dif = - Intelec[resi] + Intelec_ala[resi] - Intvdw[resi] + Intvdw_ala[resi] - IntsolvAB[resi] +\
                    IntsolvAB_ala[resi] - IntsolvA[resi] + IntsolvA_ala[resi]

        residue_ids.append(residue_id_str)
        bar_values.append(dif)
        

# Plot the bar graph
plt.bar(residue_ids, bar_values)
plt.xlabel('Residue ID')
plt.ylabel('∆∆G')
plt.title('∆∆G for Ala mutations on interface residues')
plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better readability
plt.tight_layout()  # Adjust layout to prevent clipping of labels
plt.savefig('interface_energy_changes.png')

import sys
sys.exit()
