=====
Usage
=====

To explain the usage of ``docking_py``, we will use a redocking procedure 

The same project should be launched from the ``docking`` conda environment:

.. code-block:: console

    $ conda activate docking


To use Docking Python in a project:

.. code-block:: python


    from docking_py import docking_py


Extract Ligand coordinates with ``pdb_manip_py``
------------------------------------------------


First you need to extract the ligand coordinates, we will use the `1hsg.pdb` PDB file and extract the coordinates of L-735,524 an inhibitor of the HIV proteases (resname ``MK1``) using the ``pdb_manip_py`` library (Installed with ``docking_py``):

.. code-block:: python

    >>> from pdb_manip_py import pdb_manip
    >>> # Create a Coor object
    >>> coor_1hsg = pdb_manip.Coor()
    >>> coor_1hsg.get_PDB('1hsg')
    Succeed to read file tests/input/1hsg.pdb ,  1686 atoms found
    >>> # Select res_name MK1
    >>> lig_coor = coor_1hsg.select_part_dict(
    ...     selec_dict={'res_name': 'MK1'})
    >>> # Save the ligand coordinates
    >>> lig_coor.write_pdb('lig.pdb')
    Succeed to save file lig.pdb

Extract Receptor coordinates with ``pdb_manip_py``
--------------------------------------------------

Then you need to extract the receptor coordinates, we will use the `1hsg.pdb` PDB file and extract the coordinates of the HIV II protease using the ``pdb_manip_py`` library (Installed with ``docking_py``):

.. code-block:: python

    >>> # Keep only the amino acids
    >>> rec_coor = coor_1hsg.select_part_dict(\
    ...     selec_dict={'res_name': pdb_manip.AA_DICT.keys()})
    >>> rec_coor.write_pdb('rec.pdb')
    Succeed to save file rec.pdb

Prepare Ligand and receptor structures
--------------------------------------

You need to create a ``Docking`` object, and the use the functions ``prepare_ligand()`` and ``prepare_receptor()``:

.. code-block:: python

    >>> test_dock = docking_py.Docking('test', lig_pdb='lig.pdb', rec_pdb='rec.pdb')
    >>> test_dock.prepare_ligand()
    python2.5 .../prepare_ligand4.py -l lig.pdb -B none -A\
    hydrogens -o lig.pdbqt
    >>> test_dock.prepare_receptor()
    python2.5 .../prepare_receptor4.py -r rec.pdb -A checkhydrogens\
    -o rec.pdbqt

Launch docking calculation
--------------------------

Launch the docking:

.. code-block:: python

    >>> test_dock.run_docking(name='test_dock.pdb',
                              num_modes=10,
                              energy_range=1,
                              exhaustiveness=2,
                              dock_bin='smina')

    Succeed to read file tmp/rec.pdb ,  1514 atoms found
    Grid points: [66 81 83]
    Succeed to read file tmp/rec.pdb ,  1514 atoms found
    smina --ligand lig.pdbqt --receptor rec.pdbqt --log test_dock_log.txt \
    --num_modes 10 --exhaustiveness 2 --energy_range 1 --out test_dock.pdb \
    --size_x 66.00 --size_y 81.00 --size_z 83.00 --center_x 16.07 \
    --center_y 26.49 --center_z 3.77

Analysis
--------

Extract affinity and RMSD to crystal structure:

.. code-block:: python

    >>> rmsd_list = test_dock.compute_dock_rmsd(test_dock.lig_pdbqt)
    Succeed to read file tmp/lig.pdbqt ,  50 atoms found
    Read 4 Model(s)
    Succeed to read file tmp/test_dock.pdb, 50 atoms found
    PDB file tmp/test_dock_vmd.pdb already exist, file not saved
    >>> print(rmsd_list)
    [20.730075739369596, 19.264676485215112, 14.561694090317925, 21.52437122844707]
    >>> affinity = test_dock.extract_affinity()
    >>> print(affinity)
    {1: {'affinity': -9.0, 'rmsd_low': 0.0, 'rmsd_high': 0.0}, 2: {'affinity': -8.5, 'rmsd_low': 31.972,    'rmsd_high': 35.891}, 3: {'affinity': -8.2, 'rmsd_low': 26.778, 'rmsd_high': 30.411}, 4: {'affinity': -8.1,    'rmsd_low': 31.034, 'rmsd_high': 35.639}}
