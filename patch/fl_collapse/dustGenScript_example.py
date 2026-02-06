import pyparti as dG
import numpy as np
import matplotlib.pyplot as plt

def mass_distribution_function_invGaus(r):
    return 1.1-np.exp((3*(r[0]-0.5)**2)+(3*(r[1]-0.5)**2)+(3*(r[2]-0.5)**2))
def mass_distribution_quad(r):
    return (r[0]**2+r[1]**2+r[2]**2)

ts_bins=[100,5,0.5]
n_particle=2400 
n_lagrange=2400

parameters=dG.params(n_particle,n_lagrange,mass_distribution_quad,None,ts_bins,True)
particles=dG.particle_mesh(parameters)
particles.make_setup_files()
#plt.figure()
#partx=[]
#party=[]
#for particle in particles:
#    partx.append(particle.position[0])
#    party.append(particle.position[1])
#partx=np.asarray(partx)
#party=np.asarray(party)
#plt.plot(partx,party,'.')
#plt.show()
