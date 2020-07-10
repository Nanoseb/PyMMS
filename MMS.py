#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sympy import *
from sympy.utilities.codegen import codegen
from sympy.codegen.fnodes import Module
from sympy.printing import fcode
import re

x, y, z, t = symbols('x y z t')

class MMS:
    def __init__(self, Nu=1, rho=1, 
                 U=Integer(0),
                 V=Integer(0),
                 W=Integer(0),
                 P=Integer(0),
                 Nu_t=Integer(0),
                 wall_dist=1):

        self.Nu = Nu
        self.rho = rho
        self.U = U
        self.V = V
        self.W = W
        self.P = P
        self.Nu_t = Nu_t 
        self.wall_dist = wall_dist

        self.mu = self.Nu*self.rho

        self.init_operators()
        self.init_turbulence_SA()


    def init_operators(self):     

        #  Mean strain rate Tensor
        self.Sij = Matrix([[2*diff(self.U,x)                , diff(self.U,y) + diff(self.V,x), diff(self.U,z) + diff(self.W,x)],
                           [diff(self.U,y) + diff(self.V,x) , 2*diff(self.V,y)               , diff(self.W,y) + diff(self.V,z)],
                           [diff(self.U,z) + diff(self.W,x) , diff(self.W,y) + diff(self.V,z), 2*diff(self.W,z) ]])/2


        # Material Derivative of phi
        self.MatDiff = lambda phi : diff(phi, t) + self.U*diff(phi, x) + self.V*diff(phi, y) + self.W*diff(phi, z)



    def init_turbulence_SA(self):
        """
            Initialise Spalart-Allmaras model 
        """

        c_b1 = 0.1355
        c_b2 = 0.622
        kappa = 0.41
        sigma = 2/3
        c_w2 = 0.3
        c_w3 = 2
        c_v1 = 7.1
        c_t3 = 1.2
        c_t4 = 0.5
        c_w1 = c_b1/kappa**2 + (1+c_b2)/sigma

        xi = self.Nu_t/self.Nu
        f_t2 = c_t3*exp(-c_t4*xi**2)

        self.f_v1 = xi**3/(xi**3 + c_v1**3)
        f_v2 = 1 - xi/(1+xi*self.f_v1)

        Omega = sqrt(2*(0.5*(diff(self.U, y) - diff(self.V, x))**2) + 
                     2*(0.5*(diff(self.U, z) - diff(self.W, x))**2) + 
                     2*(0.5*(diff(self.V, z) - diff(self.W, y))**2)).doit()

        S_t = Omega + self.Nu_t*f_v2/(kappa**2*self.wall_dist**2)

        #r = min(self.Nu_t/(S*kappa**2*self.wall_dist**2), 10)
        r = self.Nu_t/(S_t*kappa**2*self.wall_dist**2)
        g = (r + c_w2*(r**6 - r)).doit()

        f_w = g*((1+c_w3**6)/(g**6+c_w3**6))**(1/6)

        self.source_SA = self.MatDiff(self.Nu_t) -\
                         (c_b1*(1-f_t2)*S_t*self.Nu_t - \
                          (c_w1*f_w - c_b1*f_t2/kappa**2)*(self.Nu_t/self.wall_dist)**2 + \
                          (diff((self.Nu+self.Nu_t)*diff(self.Nu_t, x), x) + c_b2*diff(self.Nu_t, x)**2)/sigma + \
                          (diff((self.Nu+self.Nu_t)*diff(self.Nu_t, y), y) + c_b2*diff(self.Nu_t, y)**2)/sigma + \
                          (diff((self.Nu+self.Nu_t)*diff(self.Nu_t, z), z) + c_b2*diff(self.Nu_t, z)**2)/sigma)


        self.Nu_T = self.Nu_t*self.f_v1


    def get_source_SA(self):
        return self.source_SA

    def get_source_momx(self):
        return self.MatDiff(self.U) - ( -diff(self.P, x) + diff(2*self.rho*(self.Nu + self.Nu_T)*self.Sij[0,0], x) +
                                                           diff(2*self.rho*(self.Nu + self.Nu_T)*self.Sij[0,1], y) +
                                                           diff(2*self.rho*(self.Nu + self.Nu_T)*self.Sij[0,2], z)
                                     )

    def get_source_momy(self):
        return self.MatDiff(self.V) - ( -diff(self.P, y) + diff(2*self.rho*(self.Nu + self.Nu_T)*self.Sij[1,0], x) +
                                                           diff(2*self.rho*(self.Nu + self.Nu_T)*self.Sij[1,1], y) +
                                                           diff(2*self.rho*(self.Nu + self.Nu_T)*self.Sij[1,2], z)
                                     )

    def get_source_momz(self):
        return self.MatDiff(self.W) - ( -diff(self.P, z) + diff(2*self.rho*(self.Nu + self.Nu_T)*self.Sij[2,0], x) +
                                                           diff(2*self.rho*(self.Nu + self.Nu_T)*self.Sij[2,1], y) +
                                                           diff(2*self.rho*(self.Nu + self.Nu_T)*self.Sij[2,2], z)
                                     )


    def export_module(self, 
                      filename="module.F90",
                      include_field=True,
                      include_grad=True,
                      include_source=True,
                      postprocess=True):

        export_list = []

        if include_field:
            export_list.append(("U_MS", self.U.doit()))
            export_list.append(("V_MS", self.V.doit()))
            export_list.append(("W_MS", self.W.doit()))
            export_list.append(("P_MS", self.P.doit()))
            export_list.append(("Nu_t_MS", self.Nu_t.doit()))

        if include_grad:
            export_list_grad = []
            for element in export_list:
                name, equation = element
                export_list_grad.append( ("d"+name+"dx", diff(equation, x).doit()))
                export_list_grad.append( ("d"+name+"dy", diff(equation, y).doit()))
                export_list_grad.append( ("d"+name+"dz", diff(equation, z).doit()))
            export_list = export_list + export_list_grad

        if include_source:
            export_list.append(("source_MOM_U", self.get_source_momx()))
            export_list.append(("source_MOM_V", self.get_source_momy()))
            export_list.append(("source_MOM_W", self.get_source_momz()))
            export_list.append(("source_SA_Nu_t", self.get_source_SA()))


            
        [(f_name, f_code), header] = codegen(export_list, "F95", header=False, empty=True, argument_sequence=(x, y, z, t))


        module_string = fcode(Module('RANS_MMS', ['implicit none'], [f_code]), source_format='free', standard=2003)

        if postprocess:
            module_string = re.sub("\n\ *&", "", module_string)
            module_string = re.sub(" function ", " elemental function ", module_string)

        with open(filename, "w") as module_file:
            module_file.write(module_string)


        print("Function exported to {}:".format(filename))
        for element in export_list:
            print(element[0]) 






