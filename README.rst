==============
Docking Python
==============


.. image:: https://travis-ci.org/samuelmurail/docking_py.svg?branch=master
    :target: https://travis-ci.org/samuelmurail/docking_py

.. image:: https://readthedocs.org/projects/docking-py/badge/?version=latest
    :target: https://docking-py.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://codecov.io/gh/samuelmurail/docking_py/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/samuelmurail/docking_py

.. image:: https://badge.fury.io/py/docking-py.svg
    :target: https://badge.fury.io/py/docking-py

.. image:: https://static.pepy.tech/badge/docking-py
    :target: https://pepy.tech/projects/docking-py

.. image:: https://anaconda.org/bioconda/docking_py/badges/version.svg
    :target: https://anaconda.org/bioconda/docking_py

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4506970.svg
   :target: https://doi.org/10.5281/zenodo.4506970

Docking_py is a python library allowing a simplified use of the Smina, vina, qvina2 and qvinaw docking software. Docking_py can be easily automatize and scripted.


* Free software: GNU General Public License v2 (GPLv2)
* Documentation: https://docking-py.readthedocs.io.


Features
--------

* Prepare receptors and ligands.
* Launch docking using:
    * Autodock with or without **GPU** acceleration
    * Vina
    * Smina
    * Qvina2
    * Qvinaw

Citing this work
----------------

If you use the code or data in this package, please cite:

.. code-block::

    @Article{SeamDock2021,
      title = {SeamDock: An Interactive and Collaborative Online Docking Resource to Assist Small Compound Molecular Docking},
      shorttitle = {SeamDock},
      author = {Murail, Samuel and de Vries, Sjoerd J. and Rey, Julien and Moroy, Gautier and Tuff{\'e}ry, Pierre},
      year = {2021},
      journal = {Frontiers in Molecular Biosciences},
      volume = {8},
      pages = {762},
      issn = {2296-889X},
      doi = {10.3389/fmolb.2021.716466},
    }


Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
