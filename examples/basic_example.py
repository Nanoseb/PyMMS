#!/usr/bin/env python
# coding: utf-8
#
# This file is a basic example to show how to use PyMMS
# See 3D_RANS_unsteady.py for a more realistic and complex use case


from sympy import *
from PyMMS import PyMMS

# Definition of basic symbols
x, y, z, t = symbols('x y z t')

########################################
# Case initialisation
########################################
# Definition of parametric symbols
T = symbols('Period')
A = symbols('A')
B = symbols('B')

# Definition of parameters default values
global_vars = [(T, 1),
               (A, 10),
               (B, 5)]

Nu = 1
rho = 1



########################################
# Definition of field functions
########################################
U = A*cos(2*pi*t/T)*cos(x)
V = B*cos(2*pi*t/T)*cos(y)

# W verifying the continuity equation
W = - Integral(diff(U, x)+ diff(V, y), (z, 0, z))

P = Integer(1)


########################################
# MMS generation and export
########################################

# Initialisation of PyMMS with field quantities
mms = PyMMS(Nu=Nu,
            rho=rho,
            U=U,
            V=V,
            W=W,
            P=P,
            turbulence_model="none")

# Computation of source terms
mms.compute_sources()

# Export of a Fortran module file
mms.export_module("module-basic-example.F90",
                  global_vars=global_vars)


