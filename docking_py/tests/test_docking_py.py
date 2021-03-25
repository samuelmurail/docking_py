#!/usr/bin/env python

"""Tests for `docking_py` package."""

import pytest

import os
import re

from docking_py import docking
from pdb_manip_py import pdb_manip


# Test folder path
LIB_DIR = os.path.dirname(os.path.abspath(__file__))
TEST_PATH = os.path.join(LIB_DIR, "./input/")


@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string


def test_prepare_ligand_recetor(tmp_path, capsys):

    # Convert to str to avoid problem with python 3.5
    TEST_OUT = str(tmp_path)

    # Redirect pdb_manip logs
    pdb_manip.show_log()

    # Read 1hsg.pdb, extract lig.pdb and rec.pdb
    coor_1hsg = pdb_manip.Coor(os.path.join(TEST_PATH, '1hsg.pdb'))

    captured = capsys.readouterr()
    assert bool(
        re.match("Succeed to read file .+input/1hsg.pdb ,  1686 atoms found\n",
                 captured.out))

    # Extratc Protein:
    # Keep only amino acid
    rec_coor = coor_1hsg.select_part_dict(
        selec_dict={'res_name': pdb_manip.PROTEIN_RES})
    out_rec = os.path.join(TEST_OUT, 'rec.pdb')
    rec_coor.write_pdb(out_rec)

    captured = capsys.readouterr()
    assert bool(re.match("Succeed to save file .+rec.pdb\n",
                captured.out))

    # Extract ligand:
    lig_coor = coor_1hsg.select_part_dict(
        selec_dict={'res_name': 'MK1'})
    out_lig = os.path.join(TEST_OUT, 'lig.pdb')
    lig_coor.write_pdb(out_lig)

    captured = capsys.readouterr()
    assert bool(re.match("Succeed to save file .+lig.pdb\n",
                captured.out))

    # Extract ligand center of mass and maximum dimension:
    center_lig = lig_coor.center_of_mass()
    max_size = lig_coor.get_max_size()
    print("Center coordinates is {:.1f} {:.1f} {:.1f}, maximum dimension"
          " is {:.1f} Å".format(*center_lig, max_size))

    captured = capsys.readouterr()
    assert bool(re.match('Do a rotation of 99.65°\nCenter coordinates is 13.1'
                         ' 22.5 5.6, maximum dimension is 16.0 Å\n',
                captured.out))

    # Create Docking object
    test_dock = docking.Docking('test',
                                rec_pdb=out_rec,
                                lig_pdb=out_lig)

    # Prepare receptor
    test_dock.prepare_receptor()

    captured = capsys.readouterr()
    assert bool(re.match(
        ("python2.+ .+/prepare_receptor4.py -r .+/rec.pdb"
         " -A checkhydrogens -o .+/rec.pdbqt\n"),
        captured.out))

    # Check that receptor is fine
    coor_rec = pdb_manip.Coor(test_dock.rec_pdbqt)
    print('Protein atom number : {}'.format(coor_rec.select_part_dict(
        {'res_name': pdb_manip.PROTEIN_RES}).num))

    captured = capsys.readouterr()
    assert bool(re.match(
        "File name doesn't finish with .pdb read it as .pdb anyway\n"
        "Succeed to read file .+rec.pdbqt ,  1844 atoms found\n"
        "Protein atom number : 1844\n",
        captured.out))

    # Prepare Ligand
    test_dock.prepare_ligand(rigid=True)

    captured = capsys.readouterr()
    assert bool(re.match("python.+ .+prepare_ligand4.py -l lig.pdb -B "
                         "none -A hydrogens -o lig.pdbqt -Z\n",
                captured.out))

    # Check that ligand pdbqt is fine
    coor_lig = pdb_manip.Coor(test_dock.lig_pdbqt)
    print('Protein atom number : {}'.format(coor_lig.select_part_dict(
        {'res_name': pdb_manip.PROTEIN_RES}).num))

    captured = capsys.readouterr()
    assert bool(re.match(
        "File name doesn't finish with .pdb read it as .pdb anyway\n"
        "Succeed to read file .+lig.pdbqt ,  50 atoms found\n"
        "Protein atom number : 0\n",
        captured.out))


def test_smina_rigid(tmp_path, capsys):

    # Convert to str to avoid problem with python 3.5
    TEST_OUT = str(tmp_path)

    # Redirect pdb_manip logs
    pdb_manip.show_log()

    # Extract center and max_sizer:
    lig_coor = pdb_manip.Coor(os.path.join(TEST_PATH, 'lig.pdbqt'))

    center_lig = lig_coor.center_of_mass()
    max_size = lig_coor.get_max_size()
    print("Center coordinates is {:.1f} {:.1f} {:.1f}, maximum dimension"
          " is {:.1f} Å".format(*center_lig, max_size))

    captured = capsys.readouterr()
    assert bool(re.match("File name doesn't finish with .pdb read it as .pdb"
                         " anyway\n"
                         'Succeed to read file .+/input/lig.pdbqt ,  50 '
                         'atoms found\nDo a rotation of 99.71°\nCenter '
                         'coordinates is 13.1 22.5 5.5, maximum dimension is '
                         '18.0 Å\n',
                captured.out))

    # Create Docking object
    test_dock = docking.Docking(
        name='test_smina',
        rec_pdbqt=os.path.join(TEST_PATH, 'rec.pdbqt'),
        lig_pdbqt=os.path.join(TEST_PATH, 'lig.pdbqt'))

    out_dock = os.path.join(TEST_OUT, '{}_dock.pdb'.format('test_smina'))

    test_dock.run_docking(out_dock,
                          num_modes=10,
                          energy_range=10,
                          exhaustiveness=2,
                          dock_bin='smina',
                          center=center_lig,
                          grid_size=[max_size + 5] * 3,
                          seed=1)

    captured = capsys.readouterr()
    capture_line = captured.out.split('\n')
    assert bool(re.match(
        ("smina --ligand .+lig.pdbqt --receptor .+rec.pdbqt --log "
         ".+test_smina_dock_log.txt --num_modes 10 --exhaustiveness 2"
         " --energy_range 10 --out .+test_smina_dock.pdb "
         "--size_x 23.00 --size_y 23.00 --size_z 23.00"
         " --center_x 13.08 --center_y 22.52 --center_z 5.54"
         " --seed 1"),
        capture_line[0]))

    rmsd_list = test_dock.compute_dock_rmsd(test_dock.lig_pdbqt)
    assert len(rmsd_list) >= 1
    assert rmsd_list[0] < 15

    assert len(test_dock.affinity) >= 1
    assert test_dock.affinity[1]['affinity'] < -10

    # Read test_dock.dock_pdb
    coor_dock = pdb_manip.Multi_Coor(test_dock.dock_pdb)
    assert len(coor_dock.coor_list) >= 1
    assert len(coor_dock.coor_list[0].atom_dict) == 50


def test_vina_rigid(tmp_path, capsys):

    # Convert to str to avoid problem with python 3.5
    TEST_OUT = str(tmp_path)

    # Redirect pdb_manip logs
    pdb_manip.show_log()

    # Extract center and max_sizer:
    lig_coor = pdb_manip.Coor(os.path.join(TEST_PATH, 'lig.pdbqt'))

    center_lig = lig_coor.center_of_mass()
    max_size = lig_coor.get_max_size()
    print("Center coordinates is {:.1f} {:.1f} {:.1f}, maximum dimension"
          " is {:.1f} Å".format(*center_lig, max_size))

    captured = capsys.readouterr()
    assert bool(re.match("File name doesn't finish with .pdb read it as .pdb"
                         " anyway\n"
                         'Succeed to read file .+/input/lig.pdbqt ,  50 '
                         'atoms found\nDo a rotation of 99.71°\nCenter '
                         'coordinates is 13.1 22.5 5.5, maximum dimension is '
                         '18.0 Å\n',
                captured.out))

    # Create Docking object
    test_dock = docking.Docking(
        name='test_vina',
        rec_pdbqt=os.path.join(TEST_PATH, 'rec.pdbqt'),
        lig_pdbqt=os.path.join(TEST_PATH, 'lig.pdbqt'))

    out_dock = os.path.join(TEST_OUT, '{}_dock.pdb'.format('test_vina'))

    test_dock.run_docking(out_dock,
                          num_modes=10,
                          energy_range=10,
                          exhaustiveness=2,
                          dock_bin='vina',
                          center=center_lig,
                          grid_size=[max_size + 5] * 3,
                          seed=1)

    captured = capsys.readouterr()
    capture_line = captured.out.split('\n')
    assert bool(re.match(
        ("vina --ligand .+lig.pdbqt --receptor .+rec.pdbqt --log "
         ".+test_vina_dock_log.txt --num_modes 10 --exhaustiveness 2"
         " --energy_range 10 --out .+test_vina_dock.pdb "
         "--size_x 23.00 --size_y 23.00 --size_z 23.00"
         " --center_x 13.08 --center_y 22.52 --center_z 5.54"
         " --seed 1"),
        capture_line[0]))

    rmsd_list = test_dock.compute_dock_rmsd(test_dock.lig_pdbqt)
    assert len(rmsd_list) >= 1
    assert rmsd_list[0] < 15

    assert len(test_dock.affinity) >= 1
    assert test_dock.affinity[1]['affinity'] < -10

    captured = capsys.readouterr()

    test_dock.display()

    captured = capsys.readouterr()
    assert bool(re.match(
        ("name         : test_vina\n"
         "lig_pdbqt    : .+lig.pdbqt\n"
         "rec_pdbqt    : .+rec.pdbqt\n"
         "dock_pdb     : .+test_vina_dock_vmd.pdb\n"
         "dock_log     : .+test_vina_dock_log.txt\n"
         "affinity     : {1: .+'affinity': .+}}\n"),
        captured.out))

    # Read test_dock.dock_pdb
    coor_dock = pdb_manip.Multi_Coor(test_dock.dock_pdb)
    assert len(coor_dock.coor_list) >= 1
    assert len(coor_dock.coor_list[0].atom_dict) == 50


def test_autodock_rigid(tmp_path, capsys):

    # Convert to str to avoid problem with python 3.5
    TEST_OUT = str(tmp_path)

    # Redirect pdb_manip logs
    pdb_manip.show_log()

    # Extract center and max_sizer:
    lig_coor = pdb_manip.Coor(os.path.join(TEST_PATH, 'lig.pdbqt'))

    center_lig = lig_coor.center_of_mass()
    max_size = lig_coor.get_max_size()
    print("Center coordinates is {:.1f} {:.1f} {:.1f}, maximum dimension"
          " is {:.1f} Å".format(*center_lig, max_size))

    captured = capsys.readouterr()
    assert bool(re.match("File name doesn't finish with .pdb read it as .pdb"
                         " anyway\n"
                         'Succeed to read file .+/input/lig.pdbqt ,  50 '
                         'atoms found\nDo a rotation of 99.71°\nCenter '
                         'coordinates is 13.1 22.5 5.5, maximum dimension is '
                         '18.0 Å\n',
                captured.out))

    # Create Docking object
    test_dock = docking.Docking(
        name='test_autodock',
        rec_pdbqt=os.path.join(TEST_PATH, 'rec.pdbqt'),
        lig_pdbqt=os.path.join(TEST_PATH, 'lig.pdbqt'))

    # Prepare Grid
    spacing = 0.375
    test_dock.prepare_grid(out_folder=TEST_OUT, spacing=spacing,
                           center=center_lig,
                           grid_npts=[int(max_size / spacing)] * 3)

    captured = capsys.readouterr()
    assert bool(re.match(
        ("python2.+ .+prepare_gpf4.py -r rec.pdbqt -l lig.pdbqt -o "
         "test_autodock.gpf -p npts=48,48,48 -p "
         "gridcenter=13.08,22.52,5.54\nautogrid4 -p test_autodock.gpf -l "
         "test_autodock.gpf_log"),
        captured.out))

    test_dock.run_autodock(out_folder=TEST_OUT, nrun=1)

    rmsd_list = test_dock.compute_dock_rmsd(test_dock.lig_pdbqt)
    assert len(rmsd_list) >= 1
    assert rmsd_list[0] < 15
    assert test_dock.affinity[1]['affinity'] < -10

    assert bool(test_dock.lig_pdbqt.endswith("lig.pdbqt"))
    assert bool(test_dock.rec_pdbqt.endswith("rec.pdbqt"))
    assert bool(test_dock.dock_pdb.endswith("test_autodock_vmd.pdb"))
    assert bool(test_dock.gpf.endswith("test_autodock.gpf"))
    assert bool(test_dock.dock_log.endswith("test_autodock.dlg"))

    # Read test_dock.dock_pdb
    coor_dock = pdb_manip.Multi_Coor(test_dock.dock_pdb)
    assert len(coor_dock.coor_list) >= 1
    assert len(coor_dock.coor_list[0].atom_dict) == 50


def test_autodock_2_rigid(tmp_path, capsys):

    # Convert to str to avoid problem with python 3.5
    TEST_OUT = str(tmp_path)

    # Redirect pdb_manip logs
    pdb_manip.show_log()

    # Extract center and max_sizer:
    lig_coor = pdb_manip.Coor(os.path.join(TEST_PATH, 'lig.pdbqt'))

    center_lig = lig_coor.center_of_mass()
    max_size = lig_coor.get_max_size()
    print("Center coordinates is {:.1f} {:.1f} {:.1f}, maximum dimension"
          " is {:.1f} Å".format(*center_lig, max_size))

    captured = capsys.readouterr()
    assert bool(re.match("File name doesn't finish with .pdb read it as .pdb"
                         " anyway\n"
                         'Succeed to read file .+/input/lig.pdbqt ,  50 '
                         'atoms found\nDo a rotation of 99.71°\nCenter '
                         'coordinates is 13.1 22.5 5.5, maximum dimension is '
                         '18.0 Å\n',
                captured.out))

    # Create Docking object
    test_dock = docking.Docking(
        name='test_autodock_2',
        rec_pdbqt=os.path.join(TEST_PATH, 'rec.pdbqt'),
        lig_pdbqt=os.path.join(TEST_PATH, 'lig.pdbqt'))

    out_dock = os.path.join(TEST_OUT, '{}_dock.pdb'.format('test_autodock_2'))

    test_dock.run_autodock_docking(out_dock,
                                   num_modes=2,
                                   center=center_lig,
                                   grid_size=[max_size + 5] * 3)

    captured = capsys.readouterr()
    capture_line = captured.out.split('\n')
    assert bool(re.match(
        ("python2.+ .+prepare_gpf4.py -r rec.pdbqt -l lig.pdbqt -o "
         "test_autodock_2_dock.gpf -p npts=62,62,62 -p "
         "gridcenter=13.08,22.52,5.54"),
        capture_line[0]))

    assert bool(re.match(
        ("autogrid4 -p test_autodock_2_dock.gpf "
         "-l test_autodock_2_dock.gpf_log"),
        capture_line[1]))

    rmsd_list = test_dock.compute_dock_rmsd(test_dock.lig_pdbqt)
    assert len(rmsd_list) == 2
    assert rmsd_list[0] < 15

    assert len(test_dock.affinity) == 2
    assert test_dock.affinity[1]['affinity'] < -10

    assert bool(test_dock.lig_pdbqt.endswith("lig.pdbqt"))
    assert bool(test_dock.rec_pdbqt.endswith("rec.pdbqt"))
    assert bool(test_dock.dock_pdb.endswith("test_autodock_2_dock_vmd.pdb"))
    assert bool(test_dock.dock_log.endswith("test_autodock_2_dock.dlg"))
    assert bool(test_dock.gpf.endswith("test_autodock_2_dock.gpf"))

    # Read test_dock.dock_pdb
    coor_dock = pdb_manip.Multi_Coor(test_dock.dock_pdb)
    assert len(coor_dock.coor_list) == 2
    assert len(coor_dock.coor_list[0].atom_dict) == 50


def test_autodock_cpu(tmp_path, capsys):

    # Convert to str to avoid problem with python 3.5
    TEST_OUT = str(tmp_path)

    # Redirect pdb_manip logs
    pdb_manip.show_log()

    # Extract center and max_sizer:
    lig_coor = pdb_manip.Coor(os.path.join(TEST_PATH, 'lig.pdbqt'))

    center_lig = lig_coor.center_of_mass()
    max_size = lig_coor.get_max_size()
    print("Center coordinates is {:.1f} {:.1f} {:.1f}, maximum dimension"
          " is {:.1f} Å".format(*center_lig, max_size))

    captured = capsys.readouterr()
    assert bool(re.match("File name doesn't finish with .pdb read it as .pdb"
                         " anyway\n"
                         'Succeed to read file .+/input/lig.pdbqt ,  50 '
                         'atoms found\nDo a rotation of 99.71°\nCenter '
                         'coordinates is 13.1 22.5 5.5, maximum dimension is '
                         '18.0 Å\n',
                captured.out))

    # Create Docking object
    test_dock = docking.Docking(
        name='test_autodock',
        rec_pdbqt=os.path.join(TEST_PATH, 'rec.pdbqt'),
        lig_pdbqt=os.path.join(TEST_PATH, 'lig.pdbqt'))

    # Prepare Grid
    spacing = 0.375
    test_dock.prepare_grid(out_folder=TEST_OUT, spacing=spacing,
                           center=center_lig,
                           grid_npts=[int(max_size / spacing)] * 3)

    captured = capsys.readouterr()
    assert bool(re.match(
        ("python2.+ .+prepare_gpf4.py -r rec.pdbqt -l lig.pdbqt -o "
         "test_autodock.gpf -p npts=48,48,48 -p "
         "gridcenter=13.08,22.52,5.54\nautogrid4 -p test_autodock.gpf -l "
         "test_autodock.gpf_log"),
        captured.out))

    test_dock.run_autodock_cpu(out_folder=TEST_OUT, nrun=2)

    rmsd_list = test_dock.compute_dock_rmsd(test_dock.lig_pdbqt)

    assert len(rmsd_list) == 2
    assert rmsd_list[0] < 15

    assert len(test_dock.affinity) == 2
    assert test_dock.affinity[1]['affinity'] < -10

    assert bool(test_dock.lig_pdbqt.endswith("lig.pdbqt"))
    assert bool(test_dock.rec_pdbqt.endswith("rec.pdbqt"))
    assert bool(test_dock.dock_pdb.endswith("test_autodock_vmd.pdb"))
    assert bool(test_dock.dock_log.endswith("test_autodock.dlg"))
    assert bool(test_dock.gpf.endswith("test_autodock.gpf"))

    # Read test_dock.dock_pdb
    coor_dock = pdb_manip.Multi_Coor(test_dock.dock_pdb)
    assert len(coor_dock.coor_list) == 2
    assert len(coor_dock.coor_list[0].atom_dict) == 50
