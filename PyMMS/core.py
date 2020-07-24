#!/usr/bin/env python
# -*- coding: utf-8 -*-

from sympy import *
from sympy.utilities.codegen import codegen
from sympy.codegen.ast import Assignment
from sympy.codegen.fnodes import Module
from sympy.printing import fcode
import re
from .extra_models import Model_SA, Model_Menter_1eq

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
            Distance to the wall needed by "SA" turbulence model

        turbulence_model : string
            Turbulence model to choose, only "SA", "SA-noft2" and "MENTER_1eq" are implemented

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

        if self.turbulence_model == "SA-noft2":
            self.extra_models.append(Model_SA(self, noft2=True))

        if self.turbulence_model == "MENTER_1eq":
            self.extra_models.append(Model_Menter_1eq(self))

        self.source_initialized = False



    def init_operators(self):
        #  Mean strain rate Tensor
        self.Sij = Matrix([[2*diff(self.U,x)                , diff(self.U,y) + diff(self.V,x), diff(self.U,z) + diff(self.W,x)],
                           [diff(self.U,y) + diff(self.V,x) , 2*diff(self.V,y)               , diff(self.W,y) + diff(self.V,z)],
                           [diff(self.U,z) + diff(self.W,x) , diff(self.W,y) + diff(self.V,z), 2*diff(self.W,z) ]])/2


        # Material Derivative of phi
        self.MatDiff = lambda phi: diff(phi, t) + self.U*diff(phi, x) + self.V*diff(phi, y) + self.W*diff(phi, z)



    def compute_sources(self):
        """
        Computes source terms for the momentum equation and turbulence model selected. This should be run before export_module()
        """

        if self.verbose: print("Compute sources for Momentum equations")
        self.momentum_sources = [self.get_source_momx(),
                                 self.get_source_momy(),
                                 self.get_source_momz()]

        self.continuity_sources = self.get_source_continuity()

        for extra_model in self.extra_models:
            if self.verbose: print("Compute sources for {}". format(extra_model.name))
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

    def get_source_continuity(self):
        return diff(self.rho, t) + (diff(self.rho*self.U, x) +
                                    diff(self.rho*self.V, y) +
                                    diff(self.rho*self.W, z)
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
            export_list.append(("source_CONTINUITY", self.continuity_sources))

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



