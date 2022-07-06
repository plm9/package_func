.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly

    .. image:: https://api.cirrus-ci.com/github/<USER>/package_func.svg?branch=main
        :alt: Built Status
        :target: https://cirrus-ci.com/github/<USER>/package_func
    .. image:: https://readthedocs.org/projects/package_func/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://package_func.readthedocs.io/en/stable/
    .. image:: https://img.shields.io/coveralls/github/<USER>/package_func/main.svg
        :alt: Coveralls
        :target: https://coveralls.io/r/<USER>/package_func
    .. image:: https://img.shields.io/pypi/v/package_func.svg
        :alt: PyPI-Server
        :target: https://pypi.org/project/package_func/
    .. image:: https://img.shields.io/conda/vn/conda-forge/package_func.svg
        :alt: Conda-Forge
        :target: https://anaconda.org/conda-forge/package_func
    .. image:: https://pepy.tech/badge/package_func/month
        :alt: Monthly Downloads
        :target: https://pepy.tech/project/package_func
    .. image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter
        :alt: Twitter
        :target: https://twitter.com/package_func

.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/

|

============
package_func
============


    Package containing all the functions used for the probability calculations for the neutrino 
    oscillations in the vacuum. Later, we may include, functions for the oscillations in matter
    also.


All functions that are used are two files, the my_functions.py and the Prob_functions.py. 
my_functions.py contains functions such as the PMNS_param_matrix(), flavor_to_index(a,b),
D_mass(i,j) and D_mass_param(i,j). In these functions the library sympy is used and this 
wants to be changed, in order to improve the speed of the calculations.
Prob_functions.py contains all the formulas implemented for the calculations of the 
probabilities. Different approaches such as my calculations and formulas taken from 
the Review Article "Neutrino oscillations" from G.Bellini, L.Ludhova et al.


.. _pyscaffold-notes:

Note
====

This project has been set up using PyScaffold 4.2.3. For details and usage
information on PyScaffold see https://pyscaffold.org/.
