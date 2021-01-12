"""
Sensitivity analysis for a freely-propagating, premixed methane-air
flame. Computes the sensitivity of the laminar flame speed with respect
to each reaction rate constant.
"""

from __future__ import print_function

import cantera as ct
import numpy as np

# Simulation parameters
To = 300
Po = 101325
#Define the gas-mixutre and kinetics
#In this case, we are choosing a GRI3.0 gas
gas = ct.Solution('Smooke.cti')

# Create a stoichiometric CH4/Air premixed mixture 
gas.set_equivalence_ratio(1.6, 'CH4', {'O2':1.0, 'N2':3.76})
gas.TP = To, Po

width = 0.01  # m
#initial_grid = np.linspace(0,0.01 , 1000)
refine_grid = False
# Flame object
f = ct.FreeFlame(gas, width=width)
#f.FreeFlame.mdot = mdot
f.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)
f.transport_model = 'Mix'

f.solve(loglevel=1, auto=True)
print('\nmixture-averaged flamespeed = {:7f} m/s\n'.format(f.u[0]))
#print(f.multi_diff_coeffs)
#f.write_csv('CH4_adiabatic.csv', quiet=False)
#print(gas.Y)

