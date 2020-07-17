#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sympy import *
from sympy.utilities.codegen import codegen
from sympy.codegen.ast import Assignment
from sympy.codegen.fnodes import Module
from sympy.printing import fcode
import re

i = symbols('i', integer = True)
j = symbols('j', integer = True)
x, y, z, t = symbols('x y z t')


class PyMMS:
    def __init__(self, Nu=1, rho=1, 
                 U=Integer(0),
                 V=Integer(0),
                 W=Integer(0),
                 P=Integer(0),
                 Nu_t_tild=Integer(0),
                 wall_dist=1,
                 turbulence_model="SA",
                 verbose=True):

        self.Nu = Nu
        self.rho = rho
        self.U = U.doit()
        self.V = V.doit()
        self.W = W.doit()
        self.P = P.doit()
        self.Nu_t_tild = Nu_t_tild.doit()
        self.wall_dist = wall_dist.doit()
        self.verbose = verbose
        if self.verbose: print("Initialize MMS")

        self.mu = self.Nu*self.rho
        self.Nu_t = Integer(0)

        self.init_operators()

        self.extra_models = []
        self.turbulence_model = turbulence_model

        if self.turbulence_model == "SA":
            self.extra_models.append(Model_SA(self))

        if self.turbulence_model == "MENTER_1eq":
            self.extra_models.append(Model_Menter_1eq(self))

        self.source_initialized = False



    def init_operators(self):     
        #  Mean strain rate Tensor
        self.Sij = Matrix([[2*diff(self.U,x)                , diff(self.U,y) + diff(self.V,x), diff(self.U,z) + diff(self.W,x)],
                           [diff(self.U,y) + diff(self.V,x) , 2*diff(self.V,y)               , diff(self.W,y) + diff(self.V,z)],
                           [diff(self.U,z) + diff(self.W,x) , diff(self.W,y) + diff(self.V,z), 2*diff(self.W,z) ]])/2


        # Material Derivative of phi
        self.MatDiff = lambda phi : diff(phi, t) + self.U*diff(phi, x) + self.V*diff(phi, y) + self.W*diff(phi, z)



    def compute_sources(self):

        if self.verbose: print("Computing sources for Momentum equations")
        self.momentum_sources = [self.get_source_momx(),
                                 self.get_source_momy(),
                                 self.get_source_momz()]

        for extra_model in self.extra_models:
            if self.verbose: print("Computing sources for {}". format(extra_model.name))
            extra_model.compute_sources()

            # Adding extra_model source term to momentum equations
            eq_sources_momentum = extra_model.get_sources_momentum()
            self.momentum_sources = [self.momentum_sources[k] + eq_sources_momentum[k] for k in range(3)]

        self.source_initialized = True


    def get_source_momx(self):
        return self.MatDiff(self.rho*self.U) - ( -diff(self.P, x) + diff(2*self.rho*self.Nu*self.Sij[0,0], x) +
                                                                    diff(2*self.rho*self.Nu*self.Sij[0,1], y) +
                                                                    diff(2*self.rho*self.Nu*self.Sij[0,2], z)
                                                )

    def get_source_momy(self):
        return self.MatDiff(self.rho*self.V) - ( -diff(self.P, y) + diff(2*self.rho*self.Nu*self.Sij[1,0], x) +
                                                                    diff(2*self.rho*self.Nu*self.Sij[1,1], y) +
                                                                    diff(2*self.rho*self.Nu*self.Sij[1,2], z)
                                                )

    def get_source_momz(self):
        return self.MatDiff(self.rho*self.W) - ( -diff(self.P, z) + diff(2*self.rho*self.Nu*self.Sij[2,0], x) +
                                                                    diff(2*self.rho*self.Nu*self.Sij[2,1], y) +
                                                                    diff(2*self.rho*self.Nu*self.Sij[2,2], z)
                                                )


    def export_module(self, 
                      filename="module.F90",
                      global_vars=[],
                      include_field=True,
                      include_grad=True,
                      include_source=True,
                      postprocess=True,
                      suffix="_MS"):

        if self.verbose: print("Start export")
        export_list = []

        if include_field:
            if self.verbose: print(" computing field")
            export_list.append(("U" + suffix, self.U))
            export_list.append(("V" + suffix, self.V))
            export_list.append(("W" + suffix, self.W))
            export_list.append(("P" + suffix, self.P))
            for extra_model in self.extra_models:
                for name, expr in extra_model.get_variables():
                    export_list.append((name + suffix, expr))

        if include_grad:
            export_list_grad = []
            if self.verbose: print(" computing grad")
            for element in export_list:
                name, expression = element
                export_list_grad.append( ("d"+name+"dx", diff(expression, x)))
                export_list_grad.append( ("d"+name+"dy", diff(expression, y)))
                export_list_grad.append( ("d"+name+"dz", diff(expression, z)))
            export_list = export_list + export_list_grad


        if include_source and self.source_initialized:
            if self.verbose: print(" getting source")
            export_list.append(("source_MOM_U", self.momentum_sources[0]))
            export_list.append(("source_MOM_V", self.momentum_sources[1]))
            export_list.append(("source_MOM_W", self.momentum_sources[2]))

            for extra_model in self.extra_models:
                for name, term in extra_model.get_sources():
                    export_list.append(("source_{}".format(name), term))


        if self.verbose: print(" generating code")
        
        # export global variables
        global_vars_init = ""
        for var, value in global_vars:
            if var.is_integer:
                global_vars_init += "INTEGER, parameter :: {}\n".format(fcode(Assignment(var, value), source_format='free'))
            else:
                global_vars_init += "REAL*8, parameter :: {}\n".format(fcode(Assignment(var, value), source_format='free'))


        [(f_name, f_code), header] = codegen(export_list, "F95", 
                                             header=False, 
                                             empty=True, 
                                             argument_sequence=(x, y, z, t), 
                                             global_vars= [ var[0] for var in global_vars ])

        module_string = fcode(Module('RANS_MMS', ['implicit none', global_vars_init], [f_code]), source_format='free', standard=2003)


        if postprocess:
            module_string = re.sub("\n\ *&", "", module_string)
            module_string = re.sub(" function ", " elemental function ", module_string)


        if self.verbose: print(" writing to file")
        with open(filename, "w") as module_file:
            module_file.write(module_string)


        print("\nFunctions exported to {}:".format(filename))
        for element in export_list:
            print(element[0]) 




# Turbulence models definition

class Model_SA:
    def __init__(self, mms):
        self.Nu = mms.Nu
        self.rho = mms.rho
        self.U = mms.U
        self.V = mms.V
        self.W = mms.W
        self.P = mms.P
        self.Nu_t_tild = mms.Nu_t_tild
        self.wall_dist = mms.wall_dist
        self.Sij = mms.Sij
        self.MatDiff = mms.MatDiff

        self.name = "Spalart-Allmaras"

    def compute_sources(self):
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

        xi = self.Nu_t_tild/self.Nu
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

        f_w = g*((1+c_w3**6)/(g**6+c_w3**6))**(1/6)
        # f_w = 1

        self.source_SA = self.MatDiff(self.Nu_t_tild) -\
                         (c_b1*(1-f_t2)*S_t*self.Nu_t_tild - \
                         (c_w1*f_w - c_b1*f_t2/kappa**2)*(self.Nu_t_tild/self.wall_dist)**2 + \
                         ((diff((self.Nu+self.Nu_t_tild)*diff(self.Nu_t_tild, x), x) + c_b2*diff(self.Nu_t_tild, x)**2)/sigma + \
                         (diff((self.Nu+self.Nu_t_tild)*diff(self.Nu_t_tild, y), y) + c_b2*diff(self.Nu_t_tild, y)**2)/sigma + \
                         (diff((self.Nu+self.Nu_t_tild)*diff(self.Nu_t_tild, z), z) + c_b2*diff(self.Nu_t_tild, z)**2)/sigma))


        self.Nu_t = self.Nu_t_tild*f_v1

    def get_sources(self):
        return [("SA", self.source_SA)]

    def get_variables(self):
        return [("Nu_t_tild", self.Nu_t_tild), 
                ("Nu_t", self.Nu_t)]

    def get_sources_momentum(self):
        return [- (diff(2*self.rho*self.Nu_t*self.Sij[k,0], x) +
                   diff(2*self.rho*self.Nu_t*self.Sij[k,1], y) +
                   diff(2*self.rho*self.Nu_t*self.Sij[k,2], z)
                   ) for k in range(3)]



class Model_Menter_1eq:
    def __init__(self, mms):
        self.Nu = mms.Nu
        self.rho = mms.rho
        self.U = mms.U
        self.V = mms.V
        self.W = mms.W
        self.P = mms.P
        self.Nu_t_tild = mms.Nu_t_tild
        self.Sij = mms.Sij
        self.MatDiff = mms.MatDiff

        self.name = "Menter 1eq"

    def compute_sources(self):
        kappa    = 0.41
        sigma    = 1.0
        A_plus   = 13.0
        c_1      = 0.144
        c_2      = c_1/kappa**2 + 1.0/sigma
        c_3      = 7.0

        S = (sqrt(2*Sum(Sum(self.Sij[i,j]**2, (i,0,2)), (j,0,2))))

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
        return [("MENTER_1eq", self.source_MENTER_1eq)]

    def get_variables(self):
        return [("Nu_t_tild", self.Nu_t_tild), 
                ("Nu_t", self.Nu_t)]

    def get_sources_momentum(self):
        return [- (diff(2*self.rho*self.Nu_t*self.Sij[k,0], x) +
                   diff(2*self.rho*self.Nu_t*self.Sij[k,1], y) +
                   diff(2*self.rho*self.Nu_t*self.Sij[k,2], z)
                   ) for k in range(3)]



