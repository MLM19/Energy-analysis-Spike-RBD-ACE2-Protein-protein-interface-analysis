import biobb_structure_checking
import biobb_structure_checking.constants as cts
from biobb_structure_checking.structure_checking import StructureChecking
import os
import sys

base_dir_path=biobb_structure_checking.__path__[0]
args = cts.set_defaults(base_dir_path,{'notebook':True})

base_path = '/Users/msteph.padilla/Desktop/BDBI/second_year/First/Biophysics/Project/'
args['output_format'] = "pdb"
args['keep_canonical'] = False
args['input_structure_path'] = base_path + '6m0j.pdb'
args['output_structure_path'] = base_path + '6m0j_fixed.pdb'
args['output_structure_path_charges'] = base_path + '6m0j_fixed.pdbqt'
args['time_limit'] = False
args['nocache'] = False
args['copy_input'] = False
args['build_warnings'] = False
args['debug'] = False
args['verbose'] = False
args['coords_only'] = False
args['overwrite'] = False


old_stdout = sys.stdout # backup current stdout
sys.stdout = open(os.devnull, "w")
st_c = StructureChecking(base_dir_path, args)
st_c.altloc('occupancy')
st_c.ligands('All')
st_c.rem_hydrogen()
st_c.water("yes")
st_c.chiral()
st_c.backbone('--fix_atoms All --fix_chain none --add_caps none')
st_c.fixside()
st_c.getss('all')
st_c.add_hydrogen('auto')
st_c.checkall()
st_c._save_structure(args['output_structure_path'])
st_c.rem_hydrogen('yes')

sys.stdout = old_stdout # reset old stdout
os.system('check_structure -i ' + args['output_structure_path'] + ' -o ' + args['output_structure_path_charges'] + ' add_hydrogen --add_charges --add_mode auto')
