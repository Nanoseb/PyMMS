#!/usr/bin/env python
# coding: utf-8
#
# This file is a transcription of an unsteady version of the cases presented in:
# Eça, L., Hoekstra, M., & Vaz, G. (2012). Manufactured solutions for steady-flow Reynolds-averaged Navier-Stokes solvers.
#  International Journal of Computational Fluid Dynamics, 26(5), 313–332. https://doi.org/10.1080/10618562.2012.717617
#


from sympy import *
from PyMMS import PyMMS

i = symbols('i', integer = True)
x, y, z, t = symbols('x y z t')


########################################
# Case initialisation
########################################
print("Start definitions")
case = "0"

# user defined variables:
T = symbols('Period')
Re = symbols('Re')
As_1, As_2, As_3 = symbols('As_1 As_2 As_3')
Pb = symbols('Pb')
rho = symbols('Rho')

# Definition of default values
if case == "0":
    global_vars = [(Re, 10**7),
                   (Pb, 0.25),
                   (As_1, -45),
                   (As_2, 20),
                   (As_3, 10)]
if case == "A":
    global_vars = [(Re, 10**7),
                   (Pb, 0.25),
                   (As_1, -140),
                   (As_2, 40),
                   (As_3, 16)]
if case == "B":
    global_vars = [(Re, 10**7),
                   (Pb, -0.25),
                   (As_1, -80),
                   (As_2, 20),
                   (As_3, 50)]
if case == "C":
    global_vars = [(Re, 10**8),
                   (Pb, -0.1),
                   (As_1, -200),
                   (As_2, 60),
                   (As_3, 10)]
if case == "D":
    global_vars = [(Re, 10**9),
                   (Pb, 0.1),
                   (As_1, -60),
                   (As_2, 20),
                   (As_3, 5)]

global_vars.append((T, 50))
global_vars.append((rho, 1))


# Common variable for the test case
Nu = 1/Re # L=1 and V=1

# Domain start
Xmin = 0.1

# Domain end
Xmax = 1

# Domain height (Ymin=0)
Ymax = 0.4

# Eddy viscosity at outlet
Emext = 1

As = [As_1, As_2, As_3]
Al  = Array([0.0792, 0.000063, 0.005])
Bl  = Array([0.2, 0.2, 0.2])
Alf = Array([0.35, 0.4, 0.25])
Acp = Array([500, 0.25])
Aem = Array([0.4, 0.6, 10])

A1 = Al[i]*Re**(1-Bl[i])


########################################
# Definition of field functions
########################################

###############################
# time component
Ftime = 0.2 + 0.4*(1+sin(pi*(2*t/T-0.5)))

###############################
# Ums
A1 = Al[i]*Re**(1-Bl[i])
T1 = A1*y/x**Bl[i]
 
Upms = tanh(T1)
Byus = As[0]*y/exp(As[1]*y)
Bxus = 1-tanh(As[2]*(x**2-x+0.25))
Usms = Bxus*Byus

Ums = Sum(Alf[i]*Upms, (i, 0, 2)) + Usms*sin(pi*z)**2*Ftime

###############################
# Wms
Wms = Derivative(Usms, x)*sin(2*pi*z)**2/(4*pi)*Ftime

###############################
# Vms = - Integral(Derivative(Ums, x)+ Derivative(Wms, z), (y, 0, y))
Byvs = As[0]*(y+1/As[1])/As[1]/exp(As[1]*y)
Vsms = (Byvs-As[0]/As[1]**2)*Derivative(Bxus, x)
Vpms = Bl[i]*x**(Bl[i]-1)/A1*log(Upms+1) + Bl[i]*y*(Upms-1)/x
Vms = Sum(Alf[i]*Vpms, (i, 0, 2)) + Vsms*(sin(2*pi*z)*cos(2*pi*z)+sin(pi*z)**2)*Ftime

###############################
# Cp
Pcpx = x*(x*(x/3-(Xmin+Xmax)/2)+Xmin*Xmax)+1+Xmax**3/6-0.5*Xmin*Xmax**2
Pcpy = y*y*(y/3-0.5*Ymax)+1+Ymax**3/6

Pcpsx = 1.5*(x-Xmin)*pi/(Xmax-Xmin)
Pcpsy = 0.5*pi*y/Ymax

Cpms = Acp[0]*log(Pcpx)*log(Pcpy) + Acp[1]*(cos(Pcpsx)**2 * cos(Pcpsy)**2)*sin(pi*z)**2*Ftime

###############################
# Nu_t_tild
Gx = tanh(Aem[2]*(x**2-x+0.25))
Fx = (1+(Aem[0]+Aem[1]*Gx-1)*Ftime*sin(pi*z)**2)/x**0.8
Pdl = Fx*y*Re**0.2

TWMS = diff(Ums, y).subs([(y, 0)])
Yplms = sqrt(Re)*sqrt(TWMS)*y

Nu_t_tildms = (exp(-Pdl)*(0.41*Yplms-Emext)+Emext)/Re



########################################
# MMS generation and export
########################################

filename = "module-SA.F90"
print(filename)

mms = PyMMS(Nu=Nu,
            rho=rho,
            U=Ums,
            V=Vms,
            W=Wms,
            P=Cpms,
            Nu_t_tild=Nu_t_tildms,
            wall_dist=y,
            turbulence_model="SA")

mms.compute_sources()

mms.export_module(filename,
                  global_vars=global_vars)



