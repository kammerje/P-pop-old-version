"""
P_pop, a Python package to simulate exoplanet populations based on Kepler
statistics. Information about the method can be found in this paper:
http://adsabs.harvard.edu/abs/2018A%26A...609A...4K. This library is maintained
on GitHub at https://github.com/kammerje/P_pop.

===============================================================================
Please cite Kammerer & Quanz 2018 (2018A&A...609A...4K) if you use P_pop for
your research. Thank you!
===============================================================================

This repository makes use of "forecaster"
(https://github.com/chenjj2/forecaster.git, Chen & Kipping 2017,
2017ApJ...834...17C) which must be cloned into the folder called "forecaster".

Author: Jens Kammerer
Version: 3.2.0
Last edited: 08.10.19
"""


# PREAMBLE
#==============================================================================
import matplotlib.pyplot as plt
import numpy as np


# MAIN
#==============================================================================
planets = open('planets10.txt', 'r')
planets_lines = planets.readlines()

Rp = []
Porb = []
stype = []

i = 0
for line in planets_lines:
    line_temp = line.split('\t')
    if (i == 0):
        print(line_temp)
        col_Rp = np.where(np.array(line_temp) == 'Rp')[0][0]
        col_Porb = np.where(np.array(line_temp) == 'Porb')[0][0]
        col_stype = np.where(np.array(line_temp) == 'stype')[0][0]
    else:
        Rp += [float(line_temp[col_Rp])]
        Porb += [float(line_temp[col_Porb])]
        stype += [str(line_temp[col_stype])]
    i += 1
Rp = np.array(Rp)
Porb = np.array(Porb)
stype = np.array(stype)
import pdb; pdb.set_trace()
