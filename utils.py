import numpy as np

def pie_output_file(pc, outfile):
    
    fmt = '{:16}'
    with open(outfile,'w') as f:
        f.write(fmt.format('Altitude(km)'))
        f.write(fmt.format('Pressure(bar)'))
        f.write(fmt.format('Temperature(K)'))
        f.write(fmt.format('Eddy(/cm2/s)'))
        for sp in pc.dat.species_names[pc.dat.np:-2]:
            f.write(fmt.format(sp))
        f.write('\n')

        for i in range(pc.var.z.shape[0]):
            f.write(fmt.format('%.4e'%(pc.var.z[i]/1e5)))
            f.write(fmt.format('%.4e'%(pc.wrk.pressure[i]/1e6)))
            f.write(fmt.format('%.4e'%(pc.var.temperature[i])))
            f.write(fmt.format('%.4e'%(pc.var.edd[i])))
            for sp in pc.dat.species_names[pc.dat.np:-2]:
                ind = pc.dat.species_names.index(sp)
                tmp = pc.wrk.densities[ind,i]/pc.wrk.density[i]
                f.write(fmt.format('%.4e'%(tmp)))
            f.write('\n')


class PhotochemClima():
    
    def __init__(self, pc, c, T_guess = 300):
        self.pc = pc
        self.c = c
        self.T_guess = T_guess
        self.verbose = True
        self.trop_alt = -1.0
        
    def iteration(self):
        
        # Extract mixing ratios from Photochem
        P_i = np.empty(len(self.c.species_names))
        for i,sp in enumerate(self.c.species_names):
            ind = self.pc.dat.species_names.index(sp)
            mix = self.pc.wrk.densities[ind,0]/self.pc.wrk.density[0]
            P_i[i] = mix*self.pc.var.surface_pressure*1e6
        P_surf = self.pc.var.surface_pressure*1e6
        bg_gas = self.pc.dat.species_names[self.pc.dat.nsp-1]
        
        # Run climate model
        T_surf = self.c.surface_temperature_bg_gas(P_i, P_surf, bg_gas, T_guess = self.T_guess)
        self.T_guess = T_surf
        
        # Interpolate to photochem
        T = np.interp(self.pc.var.z, self.c.z, self.c.T)
        ind = np.where(self.c.T == self.c.T_trop)[0][0]
        self.trop_alt = self.c.z[ind]
        usol = self.pc.wrk.usol.copy()
        self.pc.set_temperature(T, trop_alt=self.trop_alt)
        
        # Run photochemical model
        self.pc.initialize_stepper(usol)
        tn = 0.0
        nsteps = -1
        while tn < self.pc.var.equilibrium_time:
            tn = self.pc.step()
            nsteps += 1
        self.pc.destroy_stepper()
        
        return nsteps
    
    def iterate(self, dT_surf_threshold = 0.01):
        iteration = 1
        check_convergence = False
        while True:
            nsteps = self.iteration()
            if self.verbose:
                print('Iteration = %i    nsteps = %i    T_surf = %f'%(iteration,nsteps,self.c.T_surf))
                
            if check_convergence:
                if np.abs(T_surf_prev - self.c.T_surf) < dT_surf_threshold:
                    break
                else:
                    T_surf_prev = self.c.T_surf
            else:
                T_surf_prev = self.c.T_surf
                check_convergence = True
                 
            iteration += 1