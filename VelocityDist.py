import numpy as np
import scipy as sp
from scipy.special import erf
import math 
import matplotlib.pyplot as plt

 # velocities (m/s)
v_0 = 23e4
v_e   = 24e4
v_esc = 60e4

N_0 =   np.pi**(3/2)*v_0**2*(v_0*math.erf(v_esc/v_0) - \
        2*v_esc/np.sqrt(np.pi)*np.exp(-v_esc**2/v_0**2))
  
#Truncated Maxwell-Boltzmann Distribution (Appendix B heterostructures)
a = (np.pi**(3/2)*v_0**3)/(4*v_e*N_0)

def p1(v_z):
     return (a*(math.erf((v_e-v_z)/v_0) + math.erf((v_e+v_z)/v_0) - \
      (np.pi*v_0**2/N_0)*np.exp(-v_esc**2/v_0**2)))*(v_z < v_esc - v_e) + \
     (a*(math.erf(v_esc/v_0) + math.erf((v_e-v_z)/v_0)) - \
      (np.pi*v_0**2)/(2*N_0)*(v_e + v_esc - v_z)/v_e*np.exp(-v_esc**2/v_0**2)) *\
     (v_esc - v_e < v_z < v_esc + v_e)

omega = np.linspace(0,1.6e15)
k     = np.linspace(0,6.925e9)
m_chi = 1e4

np.seterr(divide='ignore', invalid='ignore')
[K,W] = np.meshgrid(k,omega)

X = np.zeros((len(K), len(K[0])))
Y = np.zeros((len(K), len(K[0])))
contour_func = np.zeros((len(K), len(K[0])))



for i in range(len(K)):
    for j in range(len(K[0])):
        contour_func[i][j] = W[i][j]/K[i][j]*p1(W[i][j]/K[i][j] +K[i][j]/(2*m_chi))
        X[i][j] =  K[i][j]/(v_0*m_chi)
        Y[i][j] = W[i][j]/(v_0**2*m_chi) 

plt.contourf( X, Y , contour_func, 9 )
plt.colorbar()
plt.show()