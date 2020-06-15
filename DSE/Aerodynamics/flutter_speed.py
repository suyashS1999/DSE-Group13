# -*- coding: utf-8 -*-
"""

"""

import Aeroelas_analysis as aero
import sympy as sym
import numpy as np

x = sym.Symbol("x")
lmb = sym.Symbol("t")
V = 10
rho = 1.225
s= aero.half_span

mu = (aero.GJ* (s**2))/((s**2)*4*aero.EI)

a1 = aero.A_mat[0,0]*(lmb**2) + rho*aero.B_mat[0,0]*V*lmb + (rho*aero.C_mat[0,0] + x)*(V**2)

a2 = aero.A_mat[0,1]*(lmb**2) + rho*aero.B_mat[0,1]*V*lmb + (rho*aero.C_mat[0,1])*(V**2)

a3 = aero.A_mat[1,0]*(lmb**2) + rho*aero.B_mat[1,0]*V*lmb + (rho*aero.C_mat[1,0])*(V**2)

a4 = aero.A_mat[1,1]*(lmb**2) + rho*aero.B_mat[1,1]*V*lmb + (rho*aero.C_mat[1,1] + (mu*x))*(V**2)


B = np.array([[a1,a2],[a3,a4]])

B = sym.Matrix(B)
check = B.det()
check = check.expand()


b4 = check.coeff(lmb**4)
b3 = check.coeff(lmb**3)
b2 = check.coeff(lmb**2) + check.coeff(x).coeff(lmb**2)
b1 = check.coeff(lmb) + check.coeff(lmb).coeff(x)
b0 = check.coeff(lmb**0) 

result = (b4*(b1**2)) - (b1*b2*b3) + (b0*(b3**2))


sol = sym.solve(result)
