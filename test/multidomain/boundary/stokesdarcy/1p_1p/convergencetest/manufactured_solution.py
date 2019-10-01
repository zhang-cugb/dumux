import sympy as sp
import numpy as np
from sympy import *

# divergence of a matrix
def div(M):
    return np.array([sp.diff(M[0,0],x) + sp.diff(M[0,1],y) , sp.diff(M[1,0],x) + sp.diff(M[1,1],y) ])

# gradient of a vector
def grad(v):
    return np.array([[sp.diff(v[0],x),  sp.diff(v[0],y)], [sp.diff(v[1],x),  sp.diff(v[1],y)]])

# coordinates
x, y = sp.symbols("x y")

# analytical solutions for v and p in free-flow domain
vFF = np.array([(y-1.0)*(y-1.0) + x*(y-1.0) + 3.0*x - 1.0, x*(x-1.0) - 0.5*(y-1.0)*(y-1.0) - 3.0*y + 1.0])
pFF = 2.0*x + y - 1.0

# individual terms of the Navier-Stokes eq.
vvT = np.outer(vFF, vFF)
gradV = grad(vFF)
gradVT = grad(vFF).T
pI = np.array([[pFF,  0], [0,  pFF]])

# complete momentum flux and its divergence
#momentumFlux = vvT - (gradV + gradVT) +pI
momentumFlux = - (gradV + gradVT) +pI # only Stokes
divMomentumFlux = div(momentumFlux)

# solution for source term
print("divV", sp.diff(vFF[0],x) + sp.diff(vFF[1], y))
print("divMomentumFlux x:", divMomentumFlux[0])
print("divMomentumFlux y:", divMomentumFlux[1])

# analytical solution for p in Darcy domain
pD = x*(1.0-x)*(y-1.0) + pow(y-1.0, 3)/3.0 + 2.0*x + 2.0*y + 4.0

gradPdK = np.array([sp.diff(pD,x), sp.diff(pD,y)])
vD = -gradPdK
divVd = sp.diff(vD[0],x) + sp.diff(vD[1],y)

print("divVd_", simplify(divVd))
print("v x:", simplify(vD[0]))
print("v_y", simplify(vD[1]))

print("vD[1]_y=1 * n[1]", simplify(-1*vD[1].subs(y, 1)))
