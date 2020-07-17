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
    """
    Definition and operation on a Manufactured solution. Defines the Momentum equation too.
    """
    def __init__(self, Nu=1, rho=1, 
                 U=Integer(0),
                 V=Integer(0),
                 W=Integer(0),
                 P=Integer(0),
                 Nu_t_tild=Integer(0),
                 wall_dist=Integer(1),
                 turbulence_model="SA",
                 verbose=True):
        """

        Parameters
        ----------
        Nu : float or Sympy expression
            kinematic viscosity

        rho : float or Sympy expression
            fluid density

        U : Sympy expression
            x componemt of the velocity
            
        V : Sympy expression
            y componemt of the velocity

        W : Sympy expression
            z componemt of the velocity

        P : Sympy expression
            Pressure field

        Nu_t_tild : Sympy expression
            Nu_t_tild expression for use with "SA" or "MENTER_1eq" turbulence models

        wall_dist : Sympy expression
            Distance to the wall used by "SA" turbulence model

        turbulence_model : string
            Turbulence model to choose, only "SA" and "MENTER_1eq" are implemented

        verbose : Bool
            True if progress status should be outputed

        """

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
        """
        Computes source terms for the momentum equation and turbulence model selected. This should be run before export_module()
        """

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
                      suffix="_MS"):
        """
        Exports a compilable Fortran module with the field and source terms definition.

        Parameters
        ----------
        filename : string
            Name of Fortran module file to be written

        global_vars : list of tuples (Sympy symbol, default value)
            Global variables to be exported to the Module file together with their default value. Only REAL*8 and INTEGERS are supported.

        include_field : Bool
            Export field quantities

        include_grad : Bool
            Export gradient of field quantities (need include_field=True)

        include_source : Bool
            Export equations source terms

        suffix : string
            name to be appended to variable name in the module function name, used mainly to definition overloading

        """


        if self.verbose: print("Start export")
        export_list = []

        if include_field:
            if self.verbose: print(" getting field")
            export_list.append(("U" + suffix, self.U))
            export_list.append(("V" + suffix, self.V))
            export_list.append(("W" + suffix, self.W))
            export_list.append(("P" + suffix, self.P))
            for extra_model in self.extra_models:
                for name, expr in extra_model.get_variables():
                    export_list.append((name + suffix, expr))

        if include_grad:
            export_list_grad = []
            if self.verbose: print(" computing gradients")
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


        # Code generation
        [(f_name, f_code), header] = codegen(export_list, "F95", 
                                             header=False, 
                                             empty=True, 
                                             argument_sequence=(x, y, z, t), 
                                             global_vars= [ var[0] for var in global_vars ])

        module_string = fcode(Module('RANS_MMS', ['implicit none', global_vars_init], [f_code]), source_format='free', standard=2003)


        # Fix bug with empty lines in export
        module_string = re.sub("\n\ *&", "", module_string)

        # Converts function to elemental functions
        module_string = re.sub(" function ", " elemental function ", module_string)


        if self.verbose: print(" writing to file")
        with open(filename, "w") as module_file:
            module_file.write(module_string)


        print("\nFunctions exported to {}:".format(filename))
        for element in export_list:
            print(element[0]) 




# Turbulence models definition

class Model_SA:
    """
    One equation turbulence model from Spalart-Allmaras
    Spalart, P., & Allmaras, S. (1992, January 6). A one-equation turbulence model for aerodynamic flows. 30th Aerospace Sciences Meeting and Exhibit. https://doi.org/10.2514/6.1992-439
    """
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
        Computes source terms of the model
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
        self.P = mms.P
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



