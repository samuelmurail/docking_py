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
* MGLTools for `prepare_ligand4.py` and `prepare_receptor4.py`

Use `conda en create` to create it using the `.conda.yml` file. You can overide the environmnent name using the option `--name YOUR_NAME`.

.. code-block:: console

    $ conda env create -f .conda.yml

This will create an environmnet called `docking`. You will then, need to activate the environmnent:

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
