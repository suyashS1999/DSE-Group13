from numpy import*
from matplotlib import pyplot as plt
from numpy.linalg import inv

#%% ------------- Input ------------------
# I have just copied this for now, I will clean it up by importing the txt file later
#Cl = array([0.45401197, 0.59314027, 0.67983843, 0.73606313, 0.77251623, 0.80735009,
# 0.8353416,  0.86410226, 0.8862658,  0.91311351, 0.93534546, 0.96534544,
# 1.01009516, 0.93205674, 0.99882457, 1.05348972, 1.10355137, 1.14644941,
# 1.18648774, 1.22565622, 1.24058785, 1.26768158, 1.27972376, 1.28495111,
# 1.27576291, 1.22487137, 1.10590927, 1.10590927, 1.22487137, 1.27576291,
# 1.28495111, 1.27972376, 1.26768158, 1.24058785, 1.22565622, 1.18648774,
# 1.14644941, 1.10355137, 1.05348972, 0.99882457, 0.93205674, 1.01009516,
# 0.96534544, 0.93534546, 0.91311351, 0.8862658,  0.86410226, 0.8353416,
# 0.80735009, 0.77251623, 0.73606313, 0.67983843, 0.59314027, 0.45401197]);
#span_location = array([-24.35585, -23.84732, -23.33879, -22.83026, -22.32173, -21.8132,  -21.30466,
# -20.79613, -20.28759, -19.77905, -19.27052, -18.76198, -18.25344, -17.35214,
# -16.06664, -14.78112, -13.49558, -12.21003, -10.92447,  -9.6389,   -8.35332,
#  -7.06773,  -5.78213,  -4.49653,  -3.21092,  -1.9253,   -0.63968,   0.63968,
#   1.9253,    3.21092,   4.49653,   5.78213,   7.06773,   8.35332,   9.6389,
#  10.92447,  12.21003,  13.49558,  14.78112,  16.06664,  17.35214,  18.25344,
#  18.76198,  19.27052,  19.77905,  20.28759,  20.79613,  21.30466,  21.8132,
#  22.32173,  22.83026,  23.33879,  23.84732,  24.35585]);
dir = "liftdistribution.txt";
data = genfromtxt(dir);
Cl = data[1, :];
span_location = data[0, :];

#%% ------------ Functions -----------------
def monomial_basis(nodes, f):
	N = len(nodes);
	A = zeros((N, N));
	for i in range(N):
		A[:, i] = nodes**i;
	coeff = inv(A).dot(f);
	#print(coeff);
	return coeff;

def reconstruct_f(coeff, plot_nodes):
	F = 0;
	for i in range(len(coeff)):
		F += coeff[i]*plot_nodes**i;
	return F;



coeff = monomial_basis(span_location[::2], Cl[::2]);
coeff = coeff[::-1]
F = reconstruct_f(coeff, span_location);
plt.plot(span_location, Cl, "x");
plt.plot(span_location, F);
plt.grid(True);
plt.show();

#%%------------Bending Moment Equation---------------------------
from sympy import Symbol, Poly
import numpy as np

#print(coeff)
x = Symbol('x')
V = Poly(coeff, x)
#print(V)

def V(x):
    return 1.9541242997433e-30*x**26 - 1.00298256227131e-29*x**25 - 7.3443471608949e-27*x**24 + 3.28639946033419e-26*x**23 + 1.22102281205588e-23*x**22 - 4.64804343113463e-23*x**21 - 1.18220599057669e-20*x**20 + 3.70062893573879e-20*x**19 + 7.38541957812638e-18*x**18 - 1.8087943741236e-17*x**17 - 3.11669552445891e-15*x**16 + 5.51647195222334e-15*x**15 + 9.04361729764229e-13*x**14 - 9.98867726697525e-13*x**13 - 1.80266398772825e-10*x**12 + 8.51850800005809e-11*x**11 + 2.42737687552867e-8*x**10 + 2.27083965791162e-9*x**9 - 2.13542065217398e-6*x**8 - 1.16079616842799e-6*x**7 + 0.000116247428860799*x**6 + 9.35480330016276e-5*x**5 - 0.00360746497512304*x**4 - 0.0028353781451931*x**3 + 0.0553868309130224*x**2 + 0.0250522342930375*x + 1.49300235895663
#p = Polynomial()
#def V1(t):
    #return (9*t**2 - 103.72*t)
#print(V(2))
#print(span_location[29:])

Vtab = []
Vtot = 0
xtab = np.arange(0,22.5,0.01)
#print(xtab)
 

for i in range(27,len(span_location)):
    Vnum = V(xtab[i])
    Vtot = Vtot + Vnum
    Vtab.append(Vnum)
    
#print(Vtot)

#print(span_location[27:])
#print(Vtab)

plt.plot(span_location[27:],Vtab)
plt.show()
#def M(x):
 #   return x*(1.9541242997433e-30*x**26 - 1.00298256227131e-29*x**25 - 7.3443471608949e-27*x**24 + 3.28639946033419e-26*x**23 + 1.22102281205588e-23*x**22 - 4.64804343113463e-23*x**21 - 1.18220599057669e-20*x**20 + 3.70062893573879e-20*x**19 + 7.38541957812638e-18*x**18 - 1.8087943741236e-17*x**17 - 3.11669552445891e-15*x**16 + 5.51647195222334e-15*x**15 + 9.04361729764229e-13*x**14 - 9.98867726697525e-13*x**13 - 1.80266398772825e-10*x**12 + 8.51850800005809e-11*x**11 + 2.42737687552867e-8*x**10 + 2.27083965791162e-9*x**9 - 2.13542065217398e-6*x**8 - 1.16079616842799e-6*x**7 + 0.000116247428860799*x**6 + 9.35480330016276e-5*x**5 - 0.00360746497512304*x**4 - 0.0028353781451931*x**3 + 0.0553868309130224*x**2 + 0.0250522342930375*x + 1.49300235895663)
        
#from scipy.integrate import quad

Mtab = []
xtab = []
M_nosum = []
Mtot = 0

for i in range(27,len(span_location)):
    x = span_location[i]
    M = V(span_location[i])*span_location[i]
    Mtot = Mtot + M
    M_nosum.append(M)
    Mtab.append(Mtot)
    xtab.append(x)
    
#plt.plot(xtab,Mtab)

#%%------------Equilibrium equations---------------------------
P = 10
g = 10
M_nosum.append(-P*g)
xtab.append(g)

d = 49.22/2*0.55

h = 18
H_y = 0.0981
M_nosum.append(-H_y*h)
xtab.append(h)

F_br = ((Mtot-g*P-h*H_y)/d)
A_y = Vtot-P-H_y-F_br

M_nosum.append(-F_br*d)
xtab.append(d)

#%%-----------------Adding point moments to moment diagram------------------------ 
from more_itertools import sort_together

sort=sort_together([xtab, M_nosum])    #sorting points together by location

xtab_sorted = sort[0]
Mtab_sorted = sort[1]

xtab_fin =[]
Mtab_fin = []
Mtot2 = 0

for i in range(0,len(xtab_sorted)): 
    Mtot2 = Mtot2 + Mtab_sorted[i]
    Mtab_fin.append(Mtot2)
#print(Mtab_sorted)

#print(Mtab_fin)
plt.plot(xtab_sorted,Mtab_fin)
plt.show()

#%%-----------------Moment of intertia calculation------------------------ 
















#%%-----------------Calculate tensile bending stress------------------------ 
M_max = max(Mtab_fin)*10**3

M_min = min(Mtab_fin)*10**3
print(M_min)

E = 69000
I = 545283755.8622*10**(-4)*10**(-4)*10**(-4)
y = 603.15/2*10**(-3)

sigma = M_min*y/I
sigma_yield = 324*10**6


print(sigma*10**(-6))
print(sigma_yield*10**(-6))
