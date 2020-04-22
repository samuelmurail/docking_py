#!/usr/bin/env python

"""Tests for `docking_py` package."""

import pytest
import doctest

import os
import re

from docking_py import docking_py

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


def test_smina_rigid(tmp_path, capsys):

    from pdb_manip_py import pdb_manip

    # Convert to str to avoid problem with python 3.5
    TEST_OUT = str(tmp_path)

    # Read 1hsg.pdb, extract lig.pdb and rec.pdb
    coor_1hsg = pdb_manip.Coor()
    coor_1hsg.read_pdb(os.path.join(TEST_PATH, '1hsg.pdb'))

    captured = capsys.readouterr()
    assert bool(
        re.match("Succeed to read file .+input/1hsg.pdb ,  1686 atoms found\n",
                 captured.out))

    # Extratc Protein:
    # Keep only amino acid
    rec_coor = coor_1hsg.select_part_dict(
        selec_dict={'res_name': pdb_manip.AA_DICT.keys()})
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

    # Create Docking object
    test_dock = docking_py.Docking('test',
                                   rec_pdb=out_rec,
                                   lig_pdb=out_lig)

    # Prepare receptor
    test_dock.prepare_receptor()

    captured = capsys.readouterr()
    assert bool(re.match(
        ("python2.5 .+/prepare_receptor4.py -r .+/rec.pdb"
         " -A checkhydrogens -o .+/rec.pdbqt\n"),
        captured.out))

    # Check that receptor is fine
    coor_rec = pdb_manip.Coor()
    coor_rec.read_pdb(test_dock.rec_pdbqt)

    captured = capsys.readouterr()
    assert bool(re.match(
        "Succeed to read file .+rec.pdbqt ,  1844 atoms found\n",
        captured.out))

    # Prepare Ligand
    test_dock.prepare_ligand(rigid=True)

    captured = capsys.readouterr()
    assert bool(re.match("python2.5 .+prepare_ligand4.py -l .+lig.pdbqt -Z\n",
                captured.out))

    # Check that ligand pdbqt is fine
    coor_lig = pdb_manip.Coor()
    coor_lig.read_pdb(test_dock.lig_pdbqt)

    captured = capsys.readouterr()
    assert bool(re.match(
        "Succeed to read file .+lig.pdbqt ,  50 atoms found\n",
        captured.out))

    # Launch Docking
    out_dock = os.path.join(TEST_OUT, '{}_dock.pdb'.format('test'))
    test_dock.smina_docking(out_dock,
                            num_modes=10,
                            energy_range=1,
                            exhaustiveness=2)

    captured = capsys.readouterr()
    capture_line = captured.out.split('\n')
    assert bool(re.match(
        "Succeed to read file .+rec.pdb ,  1514 atoms found",
        capture_line[0]))
    assert "Grid points: [66 81 83]" == capture_line[1]
    assert bool(re.match(
        "Succeed to read file .+rec.pdb ,  1514 atoms found",
        capture_line[2]))
    assert bool(re.match(
        ("smina -l .+lig.pdbqt -r .+rec.pdbqt --log .+test_dock_log.txt"
         " --num_modes 10 --exhaustiveness 2 --energy_range 1 -o "
         ".+test_dock.pdb --size_x 66.00 --size_y 81.00 --size_z 83.00"
         " --center_x 16.07 --center_y 26.49 --center_z 3.77"),
        capture_line[3]))

    rmsd_list = test_dock.compute_dock_rmsd(test_dock.lig_pdbqt)
    assert len(rmsd_list) >= 1
    assert rmsd_list[0] < 40

    affinity = test_dock.extract_affinity()
    assert len(affinity) >= 1
    assert affinity[1]['affinity'] < 0


print("docking_py:    \t", doctest.testmod(docking_py))
