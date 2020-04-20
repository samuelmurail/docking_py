#!/usr/bin/env python3

""" Include the Docking class
"""

__author__ = "Samuel Murail"


# standard library
import os
import urllib.request

# 3rd party packages
import numpy as np
from os_command_py import os_command
from pdb_manip_py import pdb_manip


# Autorship information
__author__ = "Samuel Murail, Damien Espana, Pierre Tuffery"
__copyright__ = "Copyright 2020, RPBS"
__credits__ = ["Damien Espana", "Pierre Tuffery", "Samuel Murail"]
__license__ = "GNU General Public License v2.0"
__version__ = "1.0.1"
__maintainer__ = "Samuel Murail"
__email__ = "samuel.murail@u-paris.fr"
__status__ = "Prototype"


# Get Path of excecutables
# Check if Readthedoc is launched skip the program path searching
on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    print('Smina cannot be found')
    SMINA_LIG = ''
    SMINA_REC = ''
    SMINA_GPF = ''
    SMINA_PYTHON = ''
else:
    SMINA_LIG = os_command.which('prepare_ligand4.py')
    print("Smina ligand script is {}".format(SMINA_LIG))
    SMINA_REC = os_command.which('prepare_receptor4.py')
    print("Smina receptor script is {}".format(SMINA_REC))

    CONDA_PREFIX = os.getenv('CONDA_PREFIX')
    # With tox CONDA_PREFIX is None
    if CONDA_PREFIX is None:
        CONDA_PREFIX = '/'.join(SMINA_REC.split('/')[:-2])

    # Find a way to fix this !!
    SMINA_GPF = os.path.join(
        CONDA_PREFIX,
        'MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py')
    print("Smina grid script is {}".format(SMINA_GPF))

    # Find python 2 from conda env
    SMINA_PYTHON = os.path.join(CONDA_PREFIX, 'bin/python2.5')
    print("Python Smina is {}".format(SMINA_PYTHON))

    SMINA_BIN = os_command.which('smina')
    print("Smina executable is {}".format(SMINA_BIN))

    VINA_BIN = os_command.which('vina')
    print("Vina executable is {}".format(VINA_BIN))

    QVINA_BIN = os_command.which('qvina2')
    print("Vina executable is {}".format(VINA_BIN))

    QVINAW_BIN = os_command.which('qvinaw')
    print("Vina executable is {}".format(VINA_BIN))

    # Use : echo $CONDA_PREFIX
    # os.environ['CONDA_PREFIX']

    # import glob
    # glob.glob('/etc/r*.conf')


class Docking:
    """ The Docking class ...
    """

    def __init__(self, name, lig_pdb=None, rec_pdb=None):
        self.name = name
        self.lig_pdb = lig_pdb
        self.rec_pdb = rec_pdb

    # @property is used to get the realtive path of this variables:
    # Usefull to print command in a shorter way
    @property
    def lig_pdb(self):
        if self._lig_pdb is not None:
            return os.path.relpath(self._lig_pdb)
        return None

    @property
    def lig_pdbqt(self):
        if self._lig_pdbqt is not None:
            return os.path.relpath(self._lig_pdbqt)
        return None

    @property
    def rec_pdb(self):
        if self._rec_pdb is not None:
            return os.path.relpath(self._rec_pdb)
        return None

    @property
    def rec_pdbqt(self):
        if self._rec_pdbqt is not None:
            return os.path.relpath(self._rec_pdbqt)
        return None

    @property
    def grid(self):
        if self._grid is not None:
            return os.path.relpath(self._grid)
        return None

    @property
    def dock_pdb(self):
        if self._dock_pdb is not None:
            return os.path.relpath(self._dock_pdb)
        return None

    @property
    def dock_log(self):
        if self._dock_log is not None:
            return os.path.relpath(self._dock_log)
        return None

    # @var.setter is used to assign the full path to this variables:
    # Usefull if the working path is changed
    @lig_pdb.setter
    def lig_pdb(self, lig_pdb):
        if lig_pdb is not None:
            self._lig_pdb = os_command.full_path_and_check(lig_pdb)
        else:
            self._lig_pdb = None

    @lig_pdbqt.setter
    def lig_pdbqt(self, lig_pdbqt):
        if lig_pdbqt is not None:
            self._lig_pdbqt = os_command.full_path_and_check(lig_pdbqt)
        else:
            self._lig_pdbqt = None

    @rec_pdb.setter
    def rec_pdb(self, rec_pdb):
        if rec_pdb is not None:
            self._rec_pdb = os_command.full_path_and_check(rec_pdb)
        else:
            self._rec_pdb = None

    @rec_pdbqt.setter
    def rec_pdbqt(self, rec_pdbqt):
        if rec_pdbqt is not None:
            self._rec_pdbqt = os_command.full_path_and_check(rec_pdbqt)
        else:
            self._rec_pdbqt = None

    @grid.setter
    def grid(self, grid):
        if grid is not None:
            self._grid = os_command.full_path_and_check(grid)
        else:
            self._grid = None

    @dock_pdb.setter
    def dock_pdb(self, dock_pdb):
        if dock_pdb is not None:
            self._dock_pdb = os_command.full_path_and_check(dock_pdb)
        else:
            self._dock_pdb = None

    @dock_log.setter
    def dock_log(self, dock_log):
        if dock_log is not None:
            self._dock_log = os_command.full_path_and_check(dock_log)
        else:
            self._dock_log = None

    def prepare_ligand(self, lig_pdbqt=None, rigid=False,
                       check_file_out=True):
        """ Ligand preparation
        """

        # If lig_pdbqt is not defined use the lig_pdb name + .pdbqt
        if lig_pdbqt is None:
            lig_pdbqt = self.lig_pdb + 'qt'

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(lig_pdbqt):
            print("prepare_ligand() not launched", lig_pdbqt, "already exist")
            self.lig_pdbqt = lig_pdbqt
            return

        option = []
        if rigid:
            option.append('-Z')

        cmd_lig = os_command.Command([SMINA_PYTHON, SMINA_LIG,
                                      "-l", self.lig_pdb,
                                      "-B", 'none',
                                      "-A", 'hydrogens',
                                      "-o", lig_pdbqt] + option)
        cmd_lig.display()
        cmd_lig.run()

        self.lig_pdbqt = lig_pdbqt
        return

    def prepare_receptor(self, rec_pdbqt=None, check_file_out=True):
        """ Receptor preparation
        """

        # If lig_pdbqt is not defined use the lig_pdb name + .pdbqt
        if rec_pdbqt is None:
            rec_pdbqt = self.rec_pdb + 'qt'

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(rec_pdbqt):
            print("prepare_receptor() not launched",
                  rec_pdbqt,
                  "already exist")

            self.rec_pdbqt = rec_pdbqt
            return

        cmd_rec = os_command.Command([SMINA_PYTHON, SMINA_REC,
                                      "-r", self.rec_pdb,
                                      "-A", 'checkhydrogens',
                                      "-o", rec_pdbqt])
        cmd_rec.display()
        cmd_rec.run()

        self.rec_pdbqt = rec_pdbqt
        return

    def prepare_grid(self, grid_out=None, check_file_out=True):
        """ Grid preparation
        """

        # If grid_out is not defined use the rec_pdbqt name + .gpf
        if grid_out is None:
            grid_out = self.name + '.gpf'

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(grid_out):
            print("prepare_grid() not launched", grid_out, "already exist")
            self.grid = grid_out
            return

        cmd_grid = os_command.Command([SMINA_PYTHON, SMINA_GPF,
                                       "-r", self.rec_pdbqt,
                                       "-l", self.lig_pdbqt,
                                       "-o", grid_out])
        cmd_grid.display()
        cmd_grid.run()

        self.grid = grid_out
        return

    def get_npts(self):
        """ Extract grid size
        """
        grid_npts = []
        with open(self.grid) as file:
            for line in file:
                if 'npts' in line:
                    grid_npts.append(int(line[5:8].strip()))
                    grid_npts.append(int(line[8:11].strip()))
                    grid_npts.append(int(line[11:15].strip()))

        self.grid_npts = grid_npts
        return

    def rec_com(self):
        """ Get center of mass
        """
        rec_com = pdb_manip.Coor()
        rec_com.read_pdb(self.rec_pdb)
        self.rec_com = rec_com.center_of_mass()
        return

    def rec_grid(self, buffer_space=30):
        """ Compute grid
        """
        rec_com = pdb_manip.Coor()
        rec_com.read_pdb(self.rec_pdb)
        self.grid_npts = (np.ceil(rec_com.get_box_dim())
                          + buffer_space).astype(int)
        return

    def smina_docking(self, out_pdb=None, log=None,
                      num_modes=100, energy_range=10,
                      exhaustiveness=16, autobox=False, check_file_out=True):
        """
        Run docking
        /home/tuffery/Work/prgs/Src/Ext/SMINA/smina.static
            --center_x 56.35  --center_y -0.27 --center_z 29.72
            --size_x 30. --size_y 30. --size_z 30.
            -l peptest.pdbqt -r 2j0tA-m1.pdbqt -o test.pdb
            --num_modes 10
        """

        # If out_pdb is not defined use the rec_pdbqt name + .gpf
        if out_pdb is None:
            out_pdb = self.name + '_dock.pdb'
        if log is None:
            log = out_pdb[:-4] + '_log.txt'

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(out_pdb):
            print("smina_docking() not launched", out_pdb, "already exist")
            self.dock_pdb = out_pdb
            self.dock_log = log
            return

        option = []
        if autobox:
            # option += ['--autobox_ligand', self.lig_pdbqt]
            option += ['--autobox_add', str(10.0)]
        else:
            self.rec_grid()
            print('pnd_manip', self.grid_npts)
            self.rec_com()
            option += ["--size_x", str(self.grid_npts[0]),
                       "--size_y", str(self.grid_npts[1]),
                       "--size_z", str(self.grid_npts[2])]

            option += ["--center_x", str(self.rec_com[0]),
                       "--center_y", str(self.rec_com[1]),
                       "--center_z", str(self.rec_com[2])]

        cmd_top = os_command.Command([SMINA_BIN,
                                      "-l", self.lig_pdbqt,
                                      "-r", self.rec_pdbqt,
                                      "--log", log,
                                      "--num_modes", str(num_modes),
                                      "--exhaustiveness", str(exhaustiveness),
                                      "--energy_range", str(energy_range),
                                      "-o", out_pdb] + option)
        cmd_top.display()
        cmd_top.run()

        self.dock_pdb = out_pdb
        self.dock_log = log
        return

    def vina_docking(self, out_pdb=None, log=None, num_modes=100,
                     dock_bin='vina', energy_range=10, exhaustiveness=16,
                     check_file_out=True, cpu=None):
        """
        Run docking
        /home/tuffery/Work/prgs/Src/Ext/SMINA/smina.static
            --center_x 56.35  --center_y -0.27 --center_z 29.72
            --size_x 30. --size_y 30. --size_z 30.
            -l peptest.pdbqt -r 2j0tA-m1.pdbqt -o test.pdb
            --num_modes 10
        """

        # If out_pdb is not defined use the rec_pdbqt name + .gpf
        if out_pdb is None:
            out_pdb = self.name + '_dock.pdb'
        if log is None:
            log = out_pdb[:-4] + '_log.txt'

        if dock_bin == 'vina':
            DOCK_BIN = VINA_BIN
        elif dock_bin == 'qvina':
            DOCK_BIN = QVINA_BIN
        elif dock_bin == 'qvinaw':
            DOCK_BIN = QVINAW_BIN
        else:
            DOCK_BIN = SMINA_BIN

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(out_pdb):
            print("vina_docking() not launched", out_pdb, "already exist")
            self.dock_pdb = out_pdb
            self.dock_log = log
            return

        option = []

        self.rec_grid()
        self.rec_com()
        option += ["--size_x", str(self.grid_npts[0]),
                   "--size_y", str(self.grid_npts[1]),
                   "--size_z", str(self.grid_npts[2])]

        option += ["--center_x", str(self.rec_com[0]),
                   "--center_y", str(self.rec_com[1]),
                   "--center_z", str(self.rec_com[2])]

        if cpu is not None:
            option += ["--cpu", str(cpu)]

        cmd_top = os_command.Command([DOCK_BIN,
                                      "--ligand", self.lig_pdbqt,
                                      "--receptor", self.rec_pdbqt,
                                      "--log", log,
                                      "--num_modes", str(num_modes),
                                      "--exhaustiveness", str(exhaustiveness),
                                      "--energy_range", str(energy_range),
                                      "--out", out_pdb] + option)
        cmd_top.display()
        cmd_top.run()

        self.dock_pdb = out_pdb
        self.dock_log = log
        return

    def compute_dock_rmsd(self, selec_dict=None):

        cryst_coor = pdb_manip.Coor()
        cryst_coor.read_pdb(self.start_lig_pdb)

        dock_coor = pdb_manip.Multi_Coor()
        dock_coor.read_pdb(self.dock_pdb)
        dock_coor.write_pdb(self.dock_pdb[:-4] + '_vmd.pdb')

        rmsd = dock_coor.compute_rmsd_to(cryst_coor)

        return rmsd

    def extract_affinity(self):

        mode_read = False

        mode_info_dict = {}

        with open(self.dock_log) as dock_log:
            for line in dock_log:
                if line.startswith('-----+------------+----------+----------'):
                    mode_read = True
                    continue
                if mode_read and not line.startswith('Writing'):
                    line_split = line.split()
                    mode_info_dict[int(line_split[0])] = {
                        'affinity': float(line_split[1]),
                        'rmsd_low': float(line_split[2]),
                        'rmsd_high': float(line_split[3])}
        return mode_info_dict

    def extract_pep_rec_pdb(self, pdb_id, rec_chain,
                            lig_chain, random_rot=True):
        """ Get pdb file from the rcsb.org website.
        - Extract receptor and ligand coordinates
        - remove alternative location
        - Keep only amino acid residues
        - Save both coordinates and add it in the object
            - self.rec_pdb
            - self.lig_pdb
        """

        # Get pdb:
        out_pdb = '{}.pdb'.format(pdb_id)
        urllib.request.urlretrieve(
            'http://files.rcsb.org/download/{}.pdb'.format(pdb_id), out_pdb)

        # Treat PDB files:
        comp_coor = pdb_manip.Coor()
        comp_coor.read_pdb('{}.pdb'.format(pdb_id))
        # Keep only amino acid
        aa_comp_coor = comp_coor.select_part_dict(selec_dict={
            'res_name': pdb_manip.AA_DICT.keys()})

        # Remove alter_loc B, C, D
        alter_loc_bcd = aa_comp_coor.get_index_selection(
            {'alter_loc': ['B', 'C', 'D']})
        aa_comp_coor.del_atom_index(index_list=alter_loc_bcd)
        aa_comp_coor.change_pdb_field(change_dict={"alter_loc": ""})

        # Extract receptor pdb
        out_rec = '{}_rec.pdb'.format(pdb_id)
        rec_coor = aa_comp_coor.select_part_dict(
            selec_dict={'chain': rec_chain})
        rec_coor.write_pdb(out_rec)
        self.rec_pdb = out_rec

        # Extract ligand pdb
        out_lig = '{}_lig.pdb'.format(pdb_id)
        lig_coor = aa_comp_coor.select_part_dict(
            selec_dict={'chain': lig_chain})
        lig_coor.write_pdb(out_lig)
        self.start_lig_pdb = out_lig

        # Add random rotation
        input_lig = '{}_input_lig.pdb'.format(pdb_id)
        if random_rot:
            tau_x, tau_y, tau_z = np.random.random_sample((3,)) * 360
            lig_coor.rotation_angle(tau_x, tau_y, tau_z)
        lig_coor.write_pdb(input_lig)
        self.lig_pdb = input_lig

    def extract_align_rec_pdb(self, pdb_id, ref_pdb, rec_chain, ref_chain):
        """ Get pdb file from the rcsb.org website.
        - Extract receptor coordinates
        - remove alternative location
        - Keep only amino acid residues
        - align structure on ref
        - Save coordinates and add it in the object
            - self.rec_pdb
        """

        # Get pdb:
        out_pdb = '{}.pdb'.format(pdb_id)
        urllib.request.urlretrieve(
            'http://files.rcsb.org/download/{}.pdb'.format(pdb_id), out_pdb)
        # Treat PDB files:
        comp_coor = pdb_manip.Coor()
        comp_coor.read_pdb('{}.pdb'.format(pdb_id))
        # Keep only amino acid
        aa_comp_coor = comp_coor.select_part_dict(
            selec_dict={'res_name': pdb_manip.AA_DICT.keys()})

        # Remove alter_loc B, C, D
        alter_loc_bcd = aa_comp_coor.get_index_selection(
            {'alter_loc': ['B', 'C', 'D']})
        aa_comp_coor.del_atom_index(index_list=alter_loc_bcd)
        aa_comp_coor.change_pdb_field(change_dict={"alter_loc": ""})

        # Extract receptor pdb
        rec_coor = aa_comp_coor.select_part_dict(
            selec_dict={'chain': rec_chain})

        # Read ref_pdb
        ref_coor = pdb_manip.Coor()
        ref_coor.read_pdb(ref_pdb)
        # Keep only amino acid
        aa_ref_coor = ref_coor.select_part_dict(
            selec_dict={'res_name': pdb_manip.AA_DICT.keys()})
        # Remove alter_loc B, C, D
        alter_loc_bcd = aa_ref_coor.get_index_selection(
            {'alter_loc': ['B', 'C', 'D']})
        aa_ref_coor.del_atom_index(index_list=alter_loc_bcd)

        rec_coor.align_seq_coor_to(
            aa_ref_coor, chain_1=rec_chain, chain_2=ref_chain)

        out_rec = '{}_rec.pdb'.format(pdb_id)
        rec_coor.write_pdb(out_rec)

        self.rec_pdb = out_rec

        return

    def extract_peplig_pdb(self, pdb_id, lig_chain, random_rot=True):
        """ Get pdb file from the rcsb.org website.
        - Extract peptide ligand coordinates
        - remove alternative location
        - Keep only amino acid residues
        - Save coordinates and add it in the object
            - self.lig_pdb
        """

        # Get pdb:
        out_pdb = '{}.pdb'.format(pdb_id)
        urllib.request.urlretrieve(
            'http://files.rcsb.org/download/{}.pdb'.format(pdb_id), out_pdb)
        # Treat PDB files:
        comp_coor = pdb_manip.Coor()
        comp_coor.read_pdb('{}.pdb'.format(pdb_id))
        # Keep only amino acid
        aa_comp_coor = comp_coor.select_part_dict(
            selec_dict={'res_name': pdb_manip.AA_DICT.keys()})

        # Remove alter_loc B, C, D
        alter_loc_bcd = aa_comp_coor.get_index_selection(
            {'alter_loc': ['B', 'C', 'D']})

        aa_comp_coor.del_atom_index(index_list=alter_loc_bcd)
        aa_comp_coor.change_pdb_field(change_dict={"alter_loc": ""})

        # Extract ligand pdb
        out_lig = '{}_lig.pdb'.format(pdb_id)
        lig_coor = aa_comp_coor.select_part_dict(
            selec_dict={'chain': lig_chain})

        lig_coor.write_pdb(out_lig)
        self.start_lig_pdb = out_lig

        input_lig = '{}_input_lig.pdb'.format(pdb_id)
        if random_rot:
            tau_x, tau_y, tau_z = np.random.random_sample((3,)) * 360
            lig_coor.rotation(tau_x, tau_y, tau_z)
        lig_coor.write_pdb(input_lig)
        self.lig_pdb = input_lig


if __name__ == "__main__":

    import doctest

    print("-Test docking_py module:")

    print("docking_py:  \t", doctest.testmod())
