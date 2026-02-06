from asyncio.format_helpers import _get_function_source
from cgi import test
import numpy as np
import inspect

class params:
    def __init__(self, n_particle, n_lagrange, distribution_r, distribution_sigma=None, ts_bins=None, colocation=True):
        '''
        DustGenerator.params(self, n_particle, n_lagrange, distribution_r, distribution_sigma=None, ts_bins=None, colocation=True)
        '''
        self.n_particle=n_particle
        self.n_lagrange=n_lagrange
        self.n_tracers=self.n_lagrange-self.n_particle
        self.distribution_r=distribution_r
        self.distribution_sigma=distribution_sigma
        self.ts_bins=ts_bins
        self.colocation=colocation
        if (self.distribution_sigma is None) and (self.ts_bins is None):
            raise ValueError('neither sigma methods have been defined')
    
    def get_function_content(self):
        code_distribution_r, line_no = inspect.getsourcelines(self.distribution_r)
        code_distribution_sigma, line_no = inspect.getsourcelines(self.distribution_sigma) 
        return (code_distribution_r,code_distribution_sigma) 

    def print_params(self,destination):
        r_code,sigma_code=self.get_function_content()
        with open(destination,'w') as file:
            file.write("## Parameter file for dust particle generator ##\n")
            file.write("## This file is written in valid python synatx as to be \"eval-able\" ##\n")
            file.write("# number of particles\n")
            file.write("n_particles=",self.n_particle)
            file.write("# number of lagrangian non-backreacting particles (the rest are tracer particles)\n") 
            file.write("n_lagrange=",self.n_lagrange)
            file.write("# function definition for distribution_r\n")
            file.write(r_code)
            file.write("\n# function definition for distribution_sigma\n")
            file.write(sigma_code)
            file.write("\nts_bins=[",self.ts_bins,"]\n")
            file.write("colocation=",self.colocation,"\n")
            file.write("## End of parameter file ##")

def to_cartesian(r):
    x=r[0]*np.cos(r[1])*np.sin(r[2])
    y=r[0]*np.sin(r[1])*np.sin(r[2])
    z=r[0]*np.cos(r[2])
    return np.asarray([x,y,z])

class particle:
    def __init__(self,ts,position):
        self.ts=ts
        self.position=position

    def make_string(self):
        print_string=str(self.ts)
        for direction in self.position:
            print_string+=" "+str(direction)
        return print_string

class particle_mesh:
    def __init__(self,params):
        self.params=params
        if self.params.distribution_sigma is None:
            self.ts_bins=True
        self.assign_particle_attributes()

    def assign_particle_attributes(self):
        lagrange_particles=[]
        tracer_particles=[]
        if (self.ts_bins is True) and (self.params.colocation is True):
            nrem=self.params.n_lagrange%len(self.params.ts_bins)
            if nrem > 0:
                print("warning n_lagrange should be a multibple of the number of bins, readjusting.")
            self.params.n_lagrange-=nrem
            for n in range(self.params.n_lagrange//len(self.params.ts_bins)):
                position=self.MC_random_position()
                for tsi in self.params.ts_bins:
                    lagrange_particles.append(particle(tsi,position))
                tracer_particles.append(particle(0,position)) 
        elif self.ts_bins is True:          
            for n in range(self.params.n_lagrange):
                position=self.MC_random_position()
                ts=self.params.ts[np.floor(n/len(self.params.ts_bins))]
                lagrange_particles.append(particle(ts,position))
                tracer_particles.append(particle(0,position))
        else:
            for n in range(self.params.n_lagrange):
                position=self.MC_random_position()
                ts=self.MC_random_ts()
                lagrange_particles.append(particle(ts,position))
        self.lagrange_particles=lagrange_particles
        if False:
            tracer_particles=[]
            for n in range(self.params.n_tracers):
                    position=self.MC_random_position()
                    ts=0
                    tracer_particles=[].append(particle(ts,position))
        self.tracer_particles=tracer_particles
    
    def make_setup_files(self, printParams=False):
        if printParams is True:
            self.params.print_params("param_record.py")
        with open("particles.dat","w") as file:
            for particle in self.lagrange_particles:
                file.write(particle.make_string())
                file.write("\n")
        with open("particle_tracer.dat","w") as file_tracers:
            for particle in self.tracer_particles:
                file_tracers.write(particle.make_string())
                file_tracers.write("\n")

    def MC_random_ts(self):
       while(True):
            ts1=np.random.random(1)
            test_chance=self.params.distribution_sigma(ts1)*np.random.random(1)
            if test_chance > np.random.random(1):
                return ts1

    def MC_random_position(self):
        while(True):
            r1=np.random.random(3)
            r1[1]*=2*np.pi
            r1[0]= np.power(r1[0], 1/3)
            r1[2]= np.arccos((2*r1[2])-1)
            r1=to_cartesian(r1)
            test_chance=self.params.distribution_r(r1)*np.random.random(1)
            if test_chance > np.random.random(1):
                return r1
    
    def MC_test_safe(self):
        for i in range(100000):
            r1=np.random.random(3)
            test_chance=self.params.distribution_r(r1)*np.random.random(1)
            if test_chance > np.random.random(1):
                return r1

