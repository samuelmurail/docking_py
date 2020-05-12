#!/usr/bin/env python3

""" Include the Docking class
"""


# standard library
import os
import sys
import logging
from shutil import copy as shutil_copy

# 3rd party packages
import numpy as np
from os_command_py import os_command
from pdb_manip_py import pdb_manip

# Logging
logger = logging.getLogger(__name__)


def show_log():
    """ To use only with Doctest !!!
    Redirect logger output to sys.stdout
    """
    # Delete all handlers
    logger.handlers = []
    # Set the logger level to INFO
    logger.setLevel(logging.INFO)
    # Add sys.sdout as handler
    logger.addHandler(logging.StreamHandler(sys.stdout))
    # Show pdb_manip Logs:
    pdb_manip.show_log()


# Get Path of excecutables
# Check if Readthedoc is launched skip the program path searching
on_rtd = os.environ.get('READTHEDOCS') == 'True'
if on_rtd:
    logger.info('Smina cannot be found')
    SMINA_LIG = ''
    SMINA_REC = ''
    SMINA_GPF = ''
    SMINA_PYTHON = ''
else:

    SMINA_BIN = os_command.which('smina')
    logger.info("Smina executable is {}".format(SMINA_BIN))

    VINA_BIN = os_command.which('vina')
    logger.info("Vina executable is {}".format(VINA_BIN))

    QVINA_BIN = os_command.which('qvina2')
    logger.info("Vina executable is {}".format(VINA_BIN))

    QVINAW_BIN = os_command.which('qvinaw')
    logger.info("Vina executable is {}".format(VINA_BIN))

    # MGLTools scripts:

    PREPARE_LIG = os_command.which('prepare_ligand4.py')
    logger.info("MGLTools ligand script is {}".format(PREPARE_LIG))
    PREPARE_REC = os_command.which('prepare_receptor4.py')
    logger.info("MGLTools receptor script is {}".format(PREPARE_REC))

    # Find python 2 from conda env
    CONDA_PREFIX = os.getenv('CONDA_PREFIX')
    # With tox CONDA_PREFIX is None
    if CONDA_PREFIX is None:
        CONDA_PREFIX = '/'.join(PREPARE_REC.split('/')[:-2])

    MGLTOOL_PYTHON = os.path.join(CONDA_PREFIX, 'bin/python2.5')
    logger.info("Python Smina is {}".format(MGLTOOL_PYTHON))

    PREPARE_GPF = os.path.join(
        CONDA_PREFIX,
        'MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py')
    logger.info("MGLTools grid script is {}".format(PREPARE_GPF))

    PREPARE_DPF = os.path.join(
        CONDA_PREFIX,
        'MGLToolsPckgs/AutoDockTools/Utilities24/prepare_dpf42.py')
    logger.info("MGLTools docking prepare script is {}".format(PREPARE_DPF))

    # Autodock part
    AUTOGRID_BIN = os_command.which('autogrid4')
    logger.info("Autogrid4 executable is {}".format(AUTOGRID_BIN))

    AUTODOCK_BIN = os_command.which('autodock4')
    logger.info("Autodock4 executable is {}".format(AUTODOCK_BIN))


# Test folder path
LIB_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.join(LIB_DIR, "./tests/input/")


class Docking:
    """Docking encapsulation class.

    This class can be used to launch vina, smina, qvina and qvinaw.

    :param name: generic name of the system
    :type name: str

    :param lig_pdb: path of the ligand coordinate file (.pdb)
    :type lig_pdb: str, optional

    :param rec_pdb: path of the receptor coordinate file (.pdb)
    :type rec_pdb: str, optional

    :param lig_pdbqt: path of the ligand coordinate file (.pdbqt)
    :type lig_pdbqt: str, optional

    :param rec_pdbqt: path of the receptor coordinate file (.pdbqt)
    :type rec_pdbqt: str, optional

    :param dock_pdb: path of the docking ligand coordinate file (.pdb)
    :type dock_pdb: str, optional

    :param dock_log: path of the docking log file (.log)
    :type dock_log: str, optional

    """

    def __init__(self, name, lig_pdb=None, rec_pdb=None,
                 lig_pdbqt=None, rec_pdbqt=None):
        self.name = name
        self.lig_pdb = lig_pdb
        self.rec_pdb = rec_pdb
        self.lig_pdbqt = lig_pdbqt
        self.rec_pdbqt = rec_pdbqt
        self.dock_pdb = None
        self.dock_log = None

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
    def dpf(self):
        if self._dpf is not None:
            return os.path.relpath(self._dpf)
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

    @dpf.setter
    def dpf(self, dpf):
        if dpf is not None:
            self._dpf = os_command.full_path_and_check(dpf)
        else:
            self._dpf = None

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
                       center=False, random_rot=False,
                       check_file_out=True):
        """ Ligand preparation to `pdbqt` format using the `prepare_ligand4.py`
        command.
        Can center the ligand, could be usefull with autodock (issues when
        x,y,z > 100 Å).

        :param lig_pdbqt: output name
        :type lig_pdbqt: str, optional, default=None

        :param rigid: Flag to define if ligand is rigid
        :type rigid: bool, optional, default=False

        :param center: Flag to define if ligand have to centered
        :type center: bool, optional, default=False

        :param check_file_out: flag to check or not if file has already
            been created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        **Object requirement(s):**

            * self.lig_pdb

        **Object field(s) changed:**

            * self.lig_pdbqt

        :Example:

        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> coor_1hsg = pdb_manip.Coor(os.path.join(TEST_PATH, '1hsg.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...tests/input/1hsg.pdb ,  1686 atoms found
        >>> lig_coor = coor_1hsg.select_part_dict(\
            selec_dict={'res_name': 'MK1'})
        >>> lig_atom_num = lig_coor.num
        >>> print('Ligand has {} atoms'.format(lig_atom_num))
        Ligand has 45 atoms
        >>> out_lig = os.path.join(TEST_OUT,'lig.pdb')
        >>> lig_coor.write_pdb(out_lig) #doctest: +ELLIPSIS
        Succeed to save file .../lig.pdb
        >>> test_dock = Docking('test', lig_pdb=out_lig)
        >>> test_dock.prepare_ligand() #doctest: +ELLIPSIS
        python2.5 .../prepare_ligand4.py -l .../lig.pdb -B none -A\
 hydrogens -o .../lig.pdbqt
        >>> coor_lig = pdb_manip.Coor(test_dock.lig_pdbqt)\
        #doctest: +ELLIPSIS
        Succeed to read file .../lig.pdbqt ,  50 atoms found
        """

        # If lig_pdbqt is not defined use the lig_pdb name + .pdbqt
        if lig_pdbqt is None:
            lig_pdbqt = self.lig_pdb + 'qt'

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(lig_pdbqt):
            logger.info("prepare_ligand() not launched {} already exist".format(lig_pdbqt))
            self.ref_lig_pdb = self.lig_pdb
            if random_rot:
                self.lig_pdb = self.lig_pdb[:-4] + '_rot.pdb'
            if center:
                self.lig_pdb = self.lig_pdb[:-4] + '_center.pdb'
            self.lig_pdbqt = lig_pdbqt
            return

        option = []
        if rigid:
            option.append('-Z')

        # Define reference pdb for ligand
        self.ref_lig_pdb = self.lig_pdb

        # Add random rot:
        if random_rot:
            lig_coor = pdb_manip.Coor(self.lig_pdb)
            tau_x, tau_y, tau_z = np.random.random_sample((3,)) * 360
            lig_coor.rotation_angle(tau_x, tau_y, tau_z)
            lig_coor.write_pdb(self.lig_pdb[:-4] + '_rot.pdb')
            self.lig_pdb = self.lig_pdb[:-4] + '_rot.pdb'

        # center ligand coordinates
        if center:
            lig_coor = pdb_manip.Coor(self.lig_pdb)
            lig_com = lig_coor.center_of_mass()
            lig_coor.translate(-lig_com)
            lig_coor.write_pdb(self.lig_pdb[:-4] + '_center.pdb')
            self.lig_pdb = self.lig_pdb[:-4] + '_center.pdb'

        cmd_lig = os_command.Command([MGLTOOL_PYTHON, PREPARE_LIG,
                                      "-l", self.lig_pdb,
                                      "-B", 'none',
                                      "-A", 'hydrogens',
                                      "-o", lig_pdbqt] + option)
        cmd_lig.display()
        cmd_lig.run()

        self.lig_pdbqt = lig_pdbqt
        return

    def prepare_receptor(self, rec_pdbqt=None, check_file_out=True):
        """ Receptor preparation to `pdbqt` format using the `prepare_receptor4.py`
        command.

        :param rec_pdbqt: output name
        :type rec_pdbqt: str, optional, default=None

        :param check_file_out: flag to check or not if file has already been
            created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True

        **Object requirement(s):**

            * self.rec_pdb

        **Object field(s) changed:**

            * self.rec_pdbqt

        :Example:

        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> coor_1hsg = pdb_manip.Coor(os.path.join(TEST_PATH, '1hsg.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file .../1hsg.pdb ,  1686 atoms found
        >>> # Keep only amino acid
        >>> rec_coor = coor_1hsg.select_part_dict(\
selec_dict={'res_name': pdb_manip.PROTEIN_AA})
        >>> out_rec = os.path.join(TEST_OUT,'rec.pdb')
        >>> rec_coor.write_pdb(out_rec) #doctest: +ELLIPSIS
        Succeed to save file .../rec.pdb
        >>> rec_atom_num = rec_coor.num
        >>> print('Receptor has {} atoms'.format(rec_atom_num))
        Receptor has 1514 atoms
        >>> test_dock = Docking('test', rec_pdb=out_rec)
        >>> test_dock.prepare_receptor() #doctest: +ELLIPSIS
        python2.5 .../prepare_receptor4.py -r .../rec.pdb -A checkhydrogens\
 -o .../rec.pdbqt
        >>> coor_rec = pdb_manip.Coor(test_dock.rec_pdbqt)\
        #doctest: +ELLIPSIS
        Succeed to read file .../rec.pdbqt ,  1844 atoms found
        """

        # If lig_pdbqt is not defined use the lig_pdb name + .pdbqt
        if rec_pdbqt is None:
            rec_pdbqt = self.rec_pdb + 'qt'

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(rec_pdbqt):
            logger.info("prepare_receptor() not launched {} "
                        "already exist".format(rec_pdbqt))

            self.rec_pdbqt = rec_pdbqt
            return

        cmd_rec = os_command.Command([MGLTOOL_PYTHON, PREPARE_REC,
                                      "-r", self.rec_pdb,
                                      "-A", 'checkhydrogens',
                                      "-o", rec_pdbqt])
        cmd_rec.display()
        cmd_rec.run()

        self.rec_pdbqt = rec_pdbqt
        return

    def prepare_grid(self, out_folder, gpf_out=None,
                     spacing=0.375, grid_npts=None, center=None,
                     check_file_out=True):
        """ Grid preparation

        Launch the ``prepare_gpf4.py`` command from MGLToolsPackage.
        And ``autogrid4``.

        """

        start_dir = os.path.abspath(".")

        # Create and go in out_folder:
        # This is necessary for the autogrid creation
        os_command.create_and_go_dir(out_folder)

        if gpf_out is None:
            gpf_out = self.name + '.gpf'

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(gpf_out):
            logger.info("prepare_grid() not launched {} "
                        "already exist".format(gpf_out))
            self.gpf = gpf_out
            os.chdir(start_dir)
            return

        # Add parameters for prepare_gpf
        option_gpf = []

        if spacing != 0.375:
            option_gpf += ['-p', 'spacing={:.2f}'.format(spacing)]

        # Add grid points
        if grid_npts is None:
            self.rec_grid(spacing=spacing)
            logger.info('Grid points: {}'.format(self.grid_npts))
            grid_npts = self.grid_npts

        # Check that not point are above 255
        # Issues have been revaled with autodock gpu with more points
        clean_npts = []
        for point in grid_npts:
            if point > 255:
                logger.warning('WARNING ! Each dimension of the grid must be'
                               ' below 256. You should rise spacing !')
                clean_npts.append(255)
            else:
                clean_npts.append(point)

        option_gpf += ['-p', 'npts={:d},{:d},{:d}'.format(*clean_npts)]

        # Add grid points
        if center is None:
            self.rec_com()
            logger.info('Center: {}'.format(self.rec_com))
            center = self.rec_com
        option_gpf += ['-p', 'gridcenter={:.2f},{:.2f},{:.2f}'.format(
            *center)]

        cmd_grid = os_command.Command([MGLTOOL_PYTHON, PREPARE_GPF,
                                       "-r", self.rec_pdbqt,
                                       "-l", self.lig_pdbqt,
                                       "-o", gpf_out] + option_gpf)
        cmd_grid.display()
        cmd_grid.run()

        # The rec.pdbqt should be in the same directory as gpf_out:
        if os.path.abspath(
                os.path.dirname(self.rec_pdbqt)) != os.path.abspath("."):
            shutil_copy(self.rec_pdbqt, os.path.abspath("."))
        # The lig.pdbqt should be in the same directory as gpf_out:
        if os.path.abspath(
                os.path.dirname(self.lig_pdbqt)) != os.path.abspath("."):
            shutil_copy(self.lig_pdbqt, os.path.abspath("."))

        grid_log = gpf_out[:-4] + '.log'
        cmd_autogrid = os_command.Command([AUTOGRID_BIN,
                                           "-p", gpf_out,
                                           "-l", grid_log])
        cmd_autogrid.display()
        cmd_autogrid.run()

        self.gpf = gpf_out

        os.chdir(start_dir)

        return

    def run_autodock_cpu(self, out_folder, dock_log=None, dock_pdb=None,
                         dpf_out=None, nrun=10, param_list=[],
                         check_file_out=True):
        """
        1. Launch the ``prepare_dpf4.py`` command from MGLToolsPackage.
        2. Launch ``autodock4``
        """

        start_dir = os.path.abspath(".")

        # Define dpf name
        if dpf_out is None:
            dpf_out = self.name + '_dock_param.dpf'
        # Define dock_log name
        if dock_log is None:
            dock_log = self.name + '_dock_log.dlg'
        # Define dock_pdb name
        if dock_pdb is None:
            dock_pdb = self.name + '_dock.pdb'

        # Create and go in out_folder:
        # Run the autodock in the same directory as the dock_log file
        os_command.create_and_go_dir(out_folder)

        # Check if output files exist:
        if (check_file_out and
                os_command.check_file_and_create_path(dock_log) and
                os_command.check_file_and_create_path(dpf_out)):
            logger.info("run_autodock_cpu() not launched {} "
                        "already exist".format(dock_log))
            self.dock_log = dock_log
            self.dock_pdb = dock_pdb
            self.extract_autodock_pdb_affinity(dock_pdb)
            os.chdir(start_dir)
            return

        # Prepare docking:
        option = []
        option += ['-p', 'ga_run={:d}'.format(nrun)]

        for parameter in param_list:
            option += ['-p', parameter]

        cmd_prep = os_command.Command([MGLTOOL_PYTHON, PREPARE_DPF,
                                       "-r", self.rec_pdbqt,
                                       "-l", self.lig_pdbqt,
                                       "-o", dpf_out] + option)
        cmd_prep.display()
        cmd_prep.run()

        # Run autodock
        cmd_dock = os_command.Command([AUTODOCK_BIN,
                                       "-p", dpf_out,
                                       "-log", dock_log])
        cmd_dock.display()
        cmd_dock.run()

        self.dock_log = dock_log
        # Exract pdb form the log file:
        self.extract_autodock_pdb_affinity(dock_pdb)

        os.chdir(start_dir)
        return

    def run_autodock_gpu(self, out_folder, dock_log=None, dock_pdb=None,
                         nrun=10, check_file_out=True):
        """
        Autodock GPU arguments:

        mandatory:
        -ffile ./input/1stp/derived/1stp_protein.maps.fld
        -lfile ./input/1stp/derived/1stp_ligand.pdbqt

        opyional:
        -nrun   # LGA runs  1
        -nev    # Score evaluations (max.) per LGA run  2500000
        -ngen   # Generations (max.) per LGA run    27000
        -lsmet  Local-search method     sw (Solis-Wets)
        -lsit   # Local-search iterations (max.)    300
        -psize  Population size     150
        -mrat   Mutation rate   2 (%)
        -crat   Crossover rate  80 (%)
        -lsrat  Local-search rate   6 (%)
        -trat   Tournament (selection) rate     60 (%)
        -resnam     Name for docking output log     "docking"
        -hsym   Handle symmetry in RMSD calc.   1
        """

        # Autodock part
        AUTODOCK_GPU_BIN = os_command.which('autodock_gpu_256wi')
        logger.info("Autodock GPU executable is {}".format(AUTODOCK_GPU_BIN))

        start_dir = os.path.abspath(".")

        # Define dock_log name
        if dock_log is None:
            dock_log = self.name + '_dock_log.dlg'
        # Define dock_pdb name
        if dock_pdb is None:
            dock_pdb = self.name + '_dock.pdb'

        # Create and go in out_folder:
        # Run the autodock in the same directory as the dock_log file
        os_command.create_and_go_dir(out_folder)

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(dock_log):
            logger.info("run_autodock_gpu() not launched {} "
                        "already exist".format(dock_log))
            self.dock_log = dock_log
            self.dock_pdb = dock_pdb
            self.extract_autodock_pdb_affinity(dock_pdb)
            os.chdir(start_dir)
            return

        self.get_gridfld()

        cmd_dock = os_command.Command([AUTODOCK_GPU_BIN,
                                       "-ffile", self.gridfld,
                                       "-lfile", self.lig_pdbqt,
                                       "-nrun", str(nrun),
                                       "-resnam", dock_log[:-4]])
        cmd_dock.display()
        cmd_dock.run()

        self.dock_log = dock_log
        # Exract pdb form the log file:
        self.extract_autodock_pdb_affinity(dock_pdb)

        os.chdir(start_dir)
        return

    def run_autodock(self, out_folder, dock_log=None,
                     nrun=10, check_file_out=True):
        """
        Run autodock with cpu or gpu if available

        """

        try:
            AUTODOCK_GPU_BIN = os_command.which('autodock_gpu_256wi')
            logger.info("Autodock GPU executable is {}".format(AUTODOCK_GPU_BIN))
            logger.info("Run Autodock GPU:")
            self.run_autodock_gpu(out_folder=out_folder, dock_log=dock_log,
                                  nrun=nrun, check_file_out=check_file_out)
        except IOError:
            logger.info("Run Autodock CPU:")
            self.run_autodock_cpu(out_folder=out_folder, dock_log=dock_log,
                                  nrun=nrun, check_file_out=check_file_out)

    def extract_autodock_pdb_affinity(self, out_pdb):
        """
        Extract pdb models from the the autodock log files.
        """
        filout = open(out_pdb, 'w')

        mode_info_dict = {}

        with open(self.dock_log) as pdbfile:
            for line in pdbfile:
                if line.startswith("DOCKED: "):
                    # print(line[8:16].strip())
                    if line[8:16].strip() in ['ATOM', 'HETATM',
                                              'MODEL', 'ENDMDL']:
                        filout.write(line[8:])
                        if line[8:16].strip() == 'MODEL':
                            model = int(line[20:])
                    if line.startswith("DOCKED: USER    Estimated Free"
                                       " Energy of Binding    ="):
                        affinity = float(line.split()[8])
                        mode_info_dict[model] = {'affinity': affinity}

        filout.write("TER\n")
        filout.close()

        self.dock_pdb = out_pdb
        self.affinity = mode_info_dict

        return

    def get_gridfld(self):
        """
        Get ``gridfld`` from the ``.gpf`` file.
        """

        with open(self.gpf) as file:
            for line in file:
                if line.startswith('gridfld'):
                    gridfld = line.split()[1]

        self.gridfld = gridfld
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
        """ Get center of mass of the receptor pdb file.
        """
        if self.rec_pdb is not None:
            rec_com = pdb_manip.Coor(self.rec_pdb)
        elif self.rec_pdbqt is not None:
            rec_com = pdb_manip.Coor(self.rec_pdbqt)
        else:
            raise IOError("No receptor file defined")

        self.rec_com = rec_com.center_of_mass()
        return self.rec_com

    def rec_grid(self, buffer_space=30, spacing=1.0):
        """ Compute grid from the receptor pdb file.
        """
        if self.rec_pdb is not None:
            rec_com = pdb_manip.Coor(self.rec_pdb)
        elif self.rec_pdbqt is not None:
            rec_com = pdb_manip.Coor(self.rec_pdbqt)
        else:
            raise IOError("No receptor file defined")

        self.grid_npts = ((np.ceil(rec_com.get_box_dim()) +
                          buffer_space) / spacing).astype(int)
        return self.grid_npts

    def run_docking(self, out_pdb, log=None, dock_bin='vina',
                    num_modes=100, energy_range=10, exhaustiveness=16,
                    cpu=None, seed=None, autobox=False,
                    center=None, grid_npts=None, min_rmsd_filter=None,
                    scoring=None, check_file_out=True):
        """
        Run docking using vina, qvina, qvinaw or smina.

        :param out_pdb: PDB output name
        :type out_pdb: str

        :param log: Log ouput name
        :type log: str, optional, default=None

        :param dock_bin: Docking software name ('vina', 'qvina', 'qvinaw',
            'smina')
        :type dock_bin: str, optional, default='vina'

        :param num_modes: maximum number of binding modes to generate
        :type num_modes: int, optional, default=100

        :param energy_range: maximum energy difference between the best binding
            mode and the worst one displayed (kcal/mol)
        :type energy_range: int, optional, default=10

        :param exhaustiveness: exhaustiveness of the global search (roughly
            proportional to time): 1+
        :type exhaustiveness: int, optional, default=16

        :param cpu: the number of CPUs to use (the default is to try
            to detect the number of CPUs or, failing that, use 1)
        :type cpu: int, optional, default=None

        :param seed: explicit random seed
        :type seed: int, optional, default=None

        :param autobox: Flag to use ligand to define the docking box
        :type autobox: bool, optional, default=False

        :param center: coordinate of the center (x, y, z, Angstroms)
        :type center: list, optional, default=None

        :param grid_npts: size in the docking box (x, y, z, Angstroms)
        :type grid_npts: list, optional, default=None

        :param check_file_out: flag to check or not if file has already been
            created. If the file is present then the command break.
        :type check_file_out: bool, optional, default=True


        **Object requirement(s):**

            * self.lig_pdbqt
            * self.rec_pdbqt

        **Object field(s) changed:**

            * self.dock_pdb
            * self.dock_log

        :Example:

        """

        # If log is not defined use out_pdb minus the '.pdb' and
        # plus '_log.txt'
        if log is None:
            log = out_pdb[:-4] + '_log.txt'

        if dock_bin == 'vina':
            DOCK_BIN = VINA_BIN
        elif dock_bin == 'qvina':
            DOCK_BIN = QVINA_BIN
        elif dock_bin == 'qvinaw':
            DOCK_BIN = QVINAW_BIN
        elif dock_bin == 'smina':
            DOCK_BIN = SMINA_BIN
        else:
            logger.error('Choose an appropriate docking software among:\n'
                         '- vina\n- qvina\n- qvinaw\n- smina\n')
            return

        # Check if output files exist:
        if check_file_out and os_command.check_file_and_create_path(out_pdb):
            logger.info("vina_docking() not launched", out_pdb, "already exist")
            self.dock_pdb = out_pdb
            self.dock_log = log
            return

        option = []
        if autobox:
            if dock_bin != 'smina':
                logger.error('autobox option is only available with smina')
                raise ValueError

            option += ['--autobox_ligand', self.lig_pdbqt]
            option += ['--autobox_add', str(10.0)]

        # Define grid size:
        if grid_npts is None:
            self.rec_grid()
            logger.info('Grid points: {}'.format(self.grid_npts))
        else:
            self.grid_npts = np.array(grid_npts).astype(int)

        option += ["--size_x", '{:.2f}'.format(self.grid_npts[0]),
                   "--size_y", '{:.2f}'.format(self.grid_npts[1]),
                   "--size_z", '{:.2f}'.format(self.grid_npts[2])]

        # Define grid center:
        if center is None:
            self.rec_com()
        else:
            self.rec_com = center

        option += ["--center_x", '{:.2f}'.format(self.rec_com[0]),
                   "--center_y", '{:.2f}'.format(self.rec_com[1]),
                   "--center_z", '{:.2f}'.format(self.rec_com[2])]

        # Define cpu number:
        if cpu is not None:
            option += ["--cpu", str(cpu)]

        # Define Seed:
        if seed is not None:
            option += ["--seed", str(seed)]

        if dock_bin == 'smina':
            if min_rmsd_filter is not None:
                option += ["--min_rmsd_filter", str(min_rmsd_filter)]
            if scoring is not None:
                option += ["--scoring", str(scoring)]

        cmd_dock = os_command.Command([DOCK_BIN,
                                       "--ligand", self.lig_pdbqt,
                                       "--receptor", self.rec_pdbqt,
                                       "--log", log,
                                       "--num_modes", str(num_modes),
                                       "--exhaustiveness", str(exhaustiveness),
                                       "--energy_range", str(energy_range),
                                       "--out", out_pdb] + option)
        cmd_dock.display()
        cmd_dock.run()

        self.dock_pdb = out_pdb
        self.dock_log = log

        return

    def compute_dock_rmsd(self, ref_lig_pdb, selec_dict={}):
        """
        Compute RMSD from docking pdb to ``ref_lig_pdb``. By
        default use all atoms for RMSD calculation.
        To use only Calpha atoms define ``selec_dict={'name':['CA']}``.

        :param ref_lig_pdb: PDB reference file
        :type ref_lig_pdb: str

        :param selec_dict: Selection for RMSD calculation
        :type selec_dict: dict, optional, default={}


        :return: RMSD list
        :rtype: list

        """

        cryst_coor = pdb_manip.Coor(ref_lig_pdb)

        dock_coor = pdb_manip.Multi_Coor(self.dock_pdb)
        dock_coor.write_pdb(self.dock_pdb[:-4] + '_vmd.pdb')

        rmsd = dock_coor.compute_rmsd_to(cryst_coor, selec_dict=selec_dict)

        return rmsd

    def extract_affinity(self):
        """
        Extract affinity from the docking ``.log`` file.

        :return: Affinity and RMSD informations as a dictionnary
        :rtype: dict

        """

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

        self.affinity = mode_info_dict
        return mode_info_dict

    def extract_lig_rec_pdb(self, coor_in, folder_out, rec_select_dict,
                            lig_select_dict):
        """
        * Extract receptor and ligand coordinates from a coor file
        * remove alternative location
        * Keep only amino acid residues
        * Save both coordinates and add it in the object


        :param pdb_id: PDB ID
        :type pdb_id: str

        :param rec_chain: Chain(s) of the receptor
        :type rec_chain: list of str

        :param lig_chain: Chain(s) of the ligand
        :type lig_chain: list of str

        **Object field(s) changed:**

            * self.rec_pdb
            * self.lig_pdb

        :Example:

        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> dock_1hsg = Docking(name='1hsg')
        >>> dock_1hsg.extract_lig_rec_pdb(os.path.join(TEST_PATH, '1hsg.pdb'),\
        TEST_OUT, {'res_name': pdb_manip.PROTEIN_AA}, {'res_name': 'MK1'})\
        #doctest: +ELLIPSIS
        Succeed to read file ...1hsg.pdb ,  1686 atoms found
        Succeed to save file ...1hsg_rec.pdb
        Succeed to save file ...1hsg_input_lig.pdb
        >>> coor_lig = pdb_manip.Coor(dock_1hsg.lig_pdb) #doctest: +ELLIPSIS
        Succeed to read file ...1hsg_input_lig.pdb ,  45 atoms found
        >>> coor_rec = pdb_manip.Coor(dock_1hsg.rec_pdb) #doctest: +ELLIPSIS
        Succeed to read file ...1hsg_rec.pdb ,  1514 atoms found

        """

        # Read pdb:
        comp_coor = pdb_manip.Coor(coor_in)

        # Remove alter_loc B, C, D
        alter_loc_bcd = comp_coor.get_index_selection(
            {'alter_loc': ['B', 'C', 'D']})
        comp_coor.del_atom_index(index_list=alter_loc_bcd)
        comp_coor.change_pdb_field(change_dict={"alter_loc": ""})

        # Extract receptor pdb
        out_rec = os.path.join(folder_out, '{}_rec.pdb'.format(self.name))
        rec_coor = comp_coor.select_part_dict(
            rec_select_dict)
        rec_coor.write_pdb(out_rec)
        self.rec_pdb = out_rec

        # Extract ligand pdb
        out_lig = os.path.join(folder_out,
                               '{}_input_lig.pdb'.format(self.name))
        lig_coor = comp_coor.select_part_dict(
            lig_select_dict)
        lig_coor.write_pdb(out_lig)
        self.lig_pdb = out_lig

    def extract_receptor(self, coor_in, folder_out,
                         rec_select_dict):
        """
        * Extract receptor coordinates
        * remove alternative location
        * Keep only amino acid residues
        * align structure on ref
        * Save coordinates and add it in the object


        :param pdb_id: PDB ID
        :type pdb_id: str

        :param ref_pdb: Reference coordinates file
        :type ref_pdb:  str

        :param rec_chain: Chain(s) of the receptor
        :type rec_chain: list of str

        :param rec_chain: Chain(s) of the reference file
        :type rec_chain: list of str

        **Object field(s) changed:**

            * self.rec_pdb

        :Example:

        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> dock_1hsg = Docking(name='1hsg')
        >>> dock_1hsg.extract_receptor(os.path.join(TEST_PATH, '1hsg.pdb'),\
        TEST_OUT, {'res_name': pdb_manip.PROTEIN_AA})\
        #doctest: +ELLIPSIS
        Succeed to read file ...1hsg.pdb ,  1686 atoms found
        Succeed to save file ...1hsg_rec.pdb
        >>> coor_rec = pdb_manip.Coor(dock_1hsg.rec_pdb) #doctest: +ELLIPSIS
        Succeed to read file ...1hsg_rec.pdb ,  1514 atoms found

        """

        # Read pdb:
        comp_coor = pdb_manip.Coor(coor_in)

        # Remove alter_loc B, C, D
        alter_loc_bcd = comp_coor.get_index_selection(
            {'alter_loc': ['B', 'C', 'D']})
        comp_coor.del_atom_index(index_list=alter_loc_bcd)
        comp_coor.change_pdb_field(change_dict={"alter_loc": ""})

        # Extract receptor pdb
        rec_coor = comp_coor.select_part_dict(rec_select_dict)

        out_rec = os.path.join(folder_out, '{}_rec.pdb'.format(self.name))
        rec_coor.write_pdb(out_rec)
        self.rec_pdb = out_rec

        return

    def align_receptor(self, ref_pdb, chain_ref=['A'], chain_rec=['A']):
        """
        Align self.rec_pdb to ref_pdb.

        :Example:

        >>> pdb_manip.show_log()
        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> dock_4yob = Docking(name='4yob')
        >>> dock_4yob.extract_receptor(os.path.join(TEST_PATH, '4yob.pdb'),\
        TEST_OUT, {'res_name': pdb_manip.PROTEIN_AA})\
        #doctest: +ELLIPSIS
        Succeed to read file ...4yob.pdb ,  916 atoms found
        Succeed to save file ...4yob_rec.pdb
        >>> dock_4yob.align_receptor(os.path.join(TEST_PATH, '1hsg.pdb'))
        Succeed to read file .../4yob_rec.pdb ,  760 atoms found
        Succeed to read file .../1hsg.pdb ,  1686 atoms found
        PQITLWKRPIVTIKIGGQLKEALLNTGADDTVFEEVNLPGRWKPKLIGGIGGFVKVRQYDQVPIEICGHKVIGTVLVGPT
        ******|**|**************|*******|**||********|*******|*******| ******\
*|*********
        PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPT
        <BLANKLINE>
        PTNVIGRNLMTQIGCTLNF
        *|*|*****|*********
        PVNIIGRNLLTQIGCTLNF
        <BLANKLINE>
        Succeed to save file ...4yob_rec.pdb
        >>> coor_holo = pdb_manip.Coor(os.path.join(TEST_PATH, '1hsg.pdb'))\
        #doctest: +ELLIPSIS
        Succeed to read file ...1hsg.pdb ,  1686 atoms found
        >>> coor_rec = pdb_manip.Coor(dock_4yob.rec_pdb) #doctest: +ELLIPSIS
        Succeed to read file ...4yob_rec.pdb ,  760 atoms found
        >>> rmsd = coor_rec.compute_rmsd_to(coor_holo,\
        selec_dict={'name': ['CA'], 'chain':['A']})
        >>> print('RMSD after alignement is {:.2f} Å'.format(rmsd))
        RMSD after alignement is 1.50 Å
        """

        # Extract receptor pdb
        rec_coor = pdb_manip.Coor(self.rec_pdb)

        # Read ref_pdb
        ref_coor = pdb_manip.Coor(ref_pdb)
        # Keep only amino acid
        aa_ref_coor = ref_coor.select_part_dict(
            selec_dict={'res_name': pdb_manip.PROTEIN_AA})
        # Remove alter_loc B, C, D
        alter_loc_bcd = aa_ref_coor.get_index_selection(
            {'alter_loc': ['B', 'C', 'D']})
        aa_ref_coor.del_atom_index(index_list=alter_loc_bcd)

        rec_coor.align_seq_coor_to(aa_ref_coor, chain_1=chain_rec,
                                   chain_2=chain_ref)
        rec_coor.write_pdb(self.rec_pdb, check_file_out=False)

        return

    def extract_ligand(self, coor_in, folder_out, lig_select_dict):
        """
        * Extract ligand coordinates
        * remove alternative location
        * Save coordinates and add it in the object

        **Object field(s) changed:**

            * self.lig_pdb

        :Example:

        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> dock_1hsg = Docking(name='1hsg')
        >>> dock_1hsg.extract_ligand(os.path.join(TEST_PATH, '1hsg.pdb'),\
        TEST_OUT, {'res_name': 'MK1'})\
        #doctest: +ELLIPSIS
        Succeed to read file ...1hsg.pdb ,  1686 atoms found
        Succeed to save file ...1hsg_lig.pdb
        >>> coor_lig = pdb_manip.Coor(dock_1hsg.lig_pdb) #doctest: +ELLIPSIS
        Succeed to read file ...1hsg_lig.pdb ,  45 atoms found

        """

        # Get pdb:
        comp_coor = pdb_manip.Coor(coor_in)

        # Remove alter_loc B, C, D
        alter_loc_bcd = comp_coor.get_index_selection(
            {'alter_loc': ['B', 'C', 'D']})

        comp_coor.del_atom_index(index_list=alter_loc_bcd)
        comp_coor.change_pdb_field(change_dict={"alter_loc": ""})

        # Extract ligand pdb
        out_lig = os.path.join(folder_out, '{}_lig.pdb'.format(self.name))
        lig_coor = comp_coor.select_part_dict(
            selec_dict=lig_select_dict)

        lig_coor.write_pdb(out_lig)
        self.lig_pdb = out_lig

    def random_rot_ligand(self):
        """
        * Do a random rotation on ligand

        **Object field(s) changed:**

            * self.lig_pdb

        :Example:

        >>> TEST_OUT = str(getfixture('tmpdir'))
        >>> dock_1hsg = Docking(name='1hsg')
        >>> dock_1hsg.extract_ligand(os.path.join(TEST_PATH, '1hsg.pdb'),\
        TEST_OUT, {'res_name': 'MK1'})\
        #doctest: +ELLIPSIS
        Succeed to read file ...1hsg.pdb ,  1686 atoms found
        Succeed to save file ...1hsg_lig.pdb
        >>> coor_lig = pdb_manip.Coor(dock_1hsg.lig_pdb) #doctest: +ELLIPSIS
        Succeed to read file ...1hsg_lig.pdb ,  45 atoms found
        >>> com_before = coor_lig.center_of_mass()
        >>> dock_1hsg.random_rot_ligand()
        Succeed to read file ...1hsg_lig.pdb ,  45 atoms found
        Succeed to save file ...1hsg_lig.pdb
        >>> coor_lig = pdb_manip.Coor(dock_1hsg.lig_pdb) #doctest: +ELLIPSIS
        Succeed to read file ...1hsg_lig.pdb ,  45 atoms found
        >>> com_after = coor_lig.center_of_mass()
        >>> print('Same center of mass after rotation :{}'.format(\
com_before==com_after))
        Same center of mass after rotation :[False False False]

        ..warning:
            The function overwrite lig_pdb coordinates.
        """

        lig_coor = pdb_manip.Coor(self.lig_pdb)

        tau_x, tau_y, tau_z = np.random.random_sample((3,)) * 360
        lig_coor.rotation_angle(tau_x, tau_y, tau_z)
        lig_coor.write_pdb(self.lig_pdb, check_file_out=False)

        return


if __name__ == "__main__":

    import doctest

    print("-Test docking_py module:")

    print("docking_py:  \t", doctest.testmod())
