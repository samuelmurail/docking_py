.. highlight:: shell

============
Installation
============

1. Create Conda Environment
---------------------------

You need to create a conda environment to be able to use:

* vina
* smina
* qvina2 and qvinaw
* MGLTools for ``prepare_ligand4.py`` and ``prepare_receptor4.py`` scripts.

Use `conda en create` to create it using the ``.conda.yml`` file. You can overide the environmnent name using the option ``--name YOUR_NAME``.

.. code-block:: console

    $ conda env create -f .conda.yml

This will create an environmnet called ``docking`` (or the name you defined). You will then, need to activate the environmnent:

.. code-block:: console

    $ conda activate docking


2. Install docking_py
---------------------

The sources for Docking Python can be downloaded from the `Github repo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/samuelmurail/docking_py

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/samuelmurail/docking_py/tarball/master

Once you have a copy of the source, you can install it with:

.. code-block:: console

    $ python setup.py install

.. _Github repo: https://github.com/samuelmurail/docking_py
.. _tarball: https://github.com/samuelmurail/docking_py/tarball/master

3. Test Installation
--------------------

To test the installation, simply use ``pytest``:

.. code-block:: console

    $ pytest
    ================================ test session starts =================================
    platform darwin -- Python 3.8.2, pytest-5.4.1, py-1.8.1, pluggy-0.13.1
    rootdir: /Users/smurail/Documents/Code/docking_py, inifile: pytest.ini
    collected 4 items
    
    docking_py/docking_py.py ..                                                    [ 50%]
    tests/test_docking_py.py ..                                                    [100%]
    
    ================================= 4 passed in 17.44s =================================

