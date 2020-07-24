#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Implementation of turbulence models used by PyMMS
# Each model is defined in its own class and should provide:
# self.name
# self.compute_sources()
# self.get_sources()
# self.get_variables()
# self.get_sources_momentum()
#

from sympy import *

i = symbols('i', integer=True)
j = symbols('j', integer=True)
x, y, z, t = symbols('x y z t')




class Model_SA:
    """
    One equation turbulence model from Spalart-Allmaras
    Spalart, P., & Allmaras, S. (1992, January 6). A one-equation turbulence model for aerodynamic flows. 30th Aerospace Sciences Meeting and Exhibit. https://doi.org/10.2514/6.1992-439

    """
    def __init__(self, mms, noft2=False):
        """
        Parameters
        ----------
        noft2 : Bool
            true if the no-ft2 version of the SA model should be used

        """
        self.Nu = mms.Nu
        self.rho = mms.rho
        self.U = mms.U
        self.V = mms.V
        self.W = mms.W
        self.Nu_t_tild = mms.Nu_t_tild
        self.wall_dist = mms.wall_dist
        self.Sij = mms.Sij
        self.MatDiff = mms.MatDiff
        self.noft2 = noft2

        self.name = "Spalart-Allmaras"


    def compute_sources(self):
        """
        Computes source terms of the model
        """

        c_b1 = 0.1355
        c_b2 = 0.622
        kappa = 0.41
        sigma = 2/3
        c_w1 = c_b1/kappa**2 + (1+c_b2)/sigma
        c_w2 = 0.3
        c_w3 = 2
        c_v1 = 7.1

        xi = self.Nu_t_tild/self.Nu

        if self.noft2:
            f_t2 = 0
        else:
            c_t3 = 1.2
            c_t4 = 0.5
            f_t2 = c_t3*exp(-c_t4*xi**2)


        f_v1 = xi**3/(xi**3 + c_v1**3)
        f_v2 = 1 - xi/(1+xi*f_v1)

        Omega = sqrt(2*(0.5*(diff(self.U, y) - diff(self.V, x))**2) +
                     2*(0.5*(diff(self.U, z) - diff(self.W, x))**2) +
                     2*(0.5*(diff(self.V, z) - diff(self.W, y))**2))

        S_t = Omega + self.Nu_t_tild*f_v2/(kappa**2*self.wall_dist**2)

        #r = min(self.Nu_t_tild/(S*kappa**2*self.wall_dist**2), 10)
        r = self.Nu_t_tild/(S_t*kappa**2*self.wall_dist**2)
        g = (r + c_w2*(r**6 - r))

        # f_w = g*((1+c_w3**6)/(g**6+c_w3**6))**(1/6)
        f_w = ((1+c_w3**6)/(1+(c_w3/g)**6))**(1/6)

        self.source_SA = self.MatDiff(self.Nu_t_tild) -\
                         (c_b1*(1-f_t2)*S_t*self.Nu_t_tild - \
                         (c_w1*f_w - c_b1*f_t2/kappa**2)*(self.Nu_t_tild/self.wall_dist)**2 + \
                         ((diff((self.Nu+self.Nu_t_tild)*diff(self.Nu_t_tild, x), x) + c_b2*diff(self.Nu_t_tild, x)**2)/sigma + \
                         (diff((self.Nu+self.Nu_t_tild)*diff(self.Nu_t_tild, y), y) + c_b2*diff(self.Nu_t_tild, y)**2)/sigma + \
                         (diff((self.Nu+self.Nu_t_tild)*diff(self.Nu_t_tild, z), z) + c_b2*diff(self.Nu_t_tild, z)**2)/sigma))

        self.Nu_t = self.Nu_t_tild*f_v1


    def get_sources(self):
        """
        Returns a list of source terms and an associated name for the export.
        """
        return [("SA", self.source_SA)]


    def get_variables(self):
        """
        Returns a list of variables used in the model and an associated name.
        """
        return [("Nu_t_tild", self.Nu_t_tild),
                ("Nu_t", self.Nu_t)]


    def get_sources_momentum(self):
        """
        Returns source terms associated with the turbulence model from the Momentum equation source terms
        """
        return [- (diff(2*self.rho*self.Nu_t*self.Sij[k,0], x) +
                   diff(2*self.rho*self.Nu_t*self.Sij[k,1], y) +
                   diff(2*self.rho*self.Nu_t*self.Sij[k,2], z)
                   ) for k in range(3)]



class Model_Menter_1eq:
    """
    Menter one equation turbulence model derived from k-ε
    Menter, F. R. (1997). Eddy Viscosity Transport Equations and Their Relation to the k-ε Model. Journal of Fluids Engineering, 119(4), 876–884. https://doi.org/10.1115/1.2819511
    """
    def __init__(self, mms):
        self.Nu = mms.Nu
        self.rho = mms.rho
        self.U = mms.U
        self.V = mms.V
        self.W = mms.W
        self.Nu_t_tild = mms.Nu_t_tild
        self.Sij = mms.Sij
        self.MatDiff = mms.MatDiff

        self.name = "Menter 1eq"


    def compute_sources(self):
        """
        Computes source terms of the model
        """
        kappa    = 0.41
        sigma    = 1.0
        A_plus   = 13.0
        c_1      = 0.144
        c_2      = c_1/kappa**2 + 1.0/sigma
        c_3      = 7.0

        S = (sqrt(2*Sum(Sum(self.Sij[i,j]**2, (i,0,2)), (j,0,2)))).simplify()

        E_ke = self.Nu_t_tild**2*(diff(S, x)**2 + diff(S, y)**2 + diff(S, z)**2)/S**2

        E_bb = diff(self.Nu_t_tild, x)**2 + \
               diff(self.Nu_t_tild, y)**2 + \
               diff(self.Nu_t_tild, z)**2

        E_1e = c_3*E_bb*tanh(E_ke/(c_3*E_bb))

        D_2 = 1 - exp(-(self.Nu_t_tild/(A_plus*kappa*self.Nu))**2)
        self.Nu_t = self.Nu_t_tild*D_2

        D_1 = (self.Nu_t + self.Nu)/(self.Nu_t_tild + self.Nu)

        self.source_MENTER_1eq = self.MatDiff(self.Nu_t_tild) - \
                                 (c_1*D_1*self.Nu_t_tild*S - c_2*E_1e +
                                  diff((self.Nu+self.Nu_t_tild/sigma)*diff(self.Nu_t_tild, x), x) +
                                  diff((self.Nu+self.Nu_t_tild/sigma)*diff(self.Nu_t_tild, y), y) +
                                  diff((self.Nu+self.Nu_t_tild/sigma)*diff(self.Nu_t_tild, z), z) )


    def get_sources(self):
        """
        Returns a list of source terms and an associated name for the export.
        """
        return [("MENTER_1eq", self.source_MENTER_1eq)]


    def get_variables(self):
        """
        Returns a list of variables used in the model and an associated name.
        """
        return [("Nu_t_tild", self.Nu_t_tild),
                ("Nu_t", self.Nu_t)]


    def get_sources_momentum(self):
        """
        returns source terms associated with the turbulence model from the Momentum equation source terms
        """
        return [- (diff(2*self.rho*self.Nu_t*self.Sij[k,0], x) +
                   diff(2*self.rho*self.Nu_t*self.Sij[k,1], y) +
                   diff(2*self.rho*self.Nu_t*self.Sij[k,2], z)
                   ) for k in range(3)]




