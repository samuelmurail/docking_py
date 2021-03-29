.. highlight:: shell

============
Installation
============

1. Get sources from the `GithubRepo`_
--------------------------------------

The sources for Docking Python can be downloaded from the `GithubRepo`_.

You can either clone the public repository:

.. code-block:: console

    $ git clone git://github.com/samuelmurail/docking_py

Or download the `tarball`_:

.. code-block:: console

    $ curl -OJL https://github.com/samuelmurail/docking_py/tarball/master

Once you have a copy of the source, switch to the ``docking_py`` directorie.

.. code-block:: console

    cd docking_py


.. _GithubRepo: https://github.com/samuelmurail/docking_py
.. _tarball: https://github.com/samuelmurail/docking_py/tarball/master


2. Create Conda Environment
---------------------------


You need to create a conda environment to be able to use:  


* vina
* smina
* qvina2 and qvinaw
* MGLTools for ``prepare_ligand4.py`` and ``prepare_receptor4.py`` scripts.
* Autodock with or without GPU *support*


Use `conda en create` to create it using the ``.conda.yml`` file. You can overide the environmnent name using the option ``--name YOUR_NAME``.

.. code-block:: console

    $ conda env create -f .conda.yml

If you use a linux OS and have a GPU card, you could try the ``autodock-gpu`` version:

.. code-block:: console

    $ conda env create -f .conda_gpu.yml


This will create an environmnet called ``docking`` or ``docking_gpu`` (or the name you defined). You will then, need to activate the environmnent:

.. code-block:: console

    $ conda activate docking


3. Install docking_py
---------------------

Once you have a copy of the source and have create a conda encironment,
you can install it with:

.. code-block:: console

    $ python setup.py install



4. Test Installation
--------------------

To test the installation, simply use ``pytest``:

.. code-block:: console

    $ pytest
    ==================================== test session starts ====================================
    platform linux -- Python 3.8.2, pytest-5.4.2, py-1.9.0, pluggy-0.13.1
    rootdir: /home/murail/Documents/Code/docking_py, inifile: pytest.ini
    plugins: cov-2.10.1
    collected 13 items

    docking_py/docking.py .......                                                         [ 53%]
    docking_py/tests/test_docking_py.py ......                                            [100%]

    ============================== 13 passed, 1 warning in 21.18s ===============================

