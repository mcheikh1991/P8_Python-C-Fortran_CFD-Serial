#!/bin/bash
python -m numpy.f2py -c oned_euler_fluxes_v5.f90 -m Euler_Flux_Schemes_1D
python -m numpy.f2py -c twod_euler_fluxes_v2.f90 -m Euler_Flux_Schemes_2D