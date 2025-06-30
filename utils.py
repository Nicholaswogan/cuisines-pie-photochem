import numpy as np
import numpy as np
import numba as nb
from numba import types
from scipy import constants as const
from scipy import integrate
from copy import deepcopy
from photochem import EvoAtmosphere, PhotoException
from tempfile import NamedTemporaryFile

def pie_output_file(pc, outfile):
    
    fmt = '{:16}'
    with open(outfile,'w') as f:
        f.write(fmt.format('Altitude(km)'))
        f.write(fmt.format('Pressure(bar)'))
        f.write(fmt.format('Temperature(K)'))
        f.write(fmt.format('Eddy(cm2/s)'))
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
                tmp = np.maximum(pc.wrk.densities[ind,i]/pc.wrk.density[i],1e-200)
                f.write(fmt.format('%.4e'%(tmp)))
            f.write('\n')

class RobustData():
    
    def __init__(self):
        self.min_mix_reset = -1e-13
        self.nsteps_total = None
        self.nerrors_total = None
        self.nerrors_before_giveup = 10

class EvoAtmosphereRobust(EvoAtmosphere):

    def __init__(self, mechanism_file, settings_file, flux_file, atmosphere_file=None, data_dir=None):

        with NamedTemporaryFile('w',suffix='.txt') as f:
            f.write(ATMOSPHERE_INIT)
            f.flush()
            name = atmosphere_file
            if name is None:
                name = f.name
            super().__init__(
                mechanism_file, 
                settings_file, 
                flux_file,
                name,
                data_dir
            )

        self.rdat = RobustData()

        # Values in photochem to adjust
        self.var.verbose = 1
        self.var.upwind_molec_diff = True
        self.var.autodiff = True
        self.var.atol = 1.0e-23
        self.var.equilibrium_time = 1e17
        self.set_particle_parameters(1, 100, 10)

    def set_particle_parameters(self, smooth_factor, k_cond, k_evap):
        for i in range(len(self.var.cond_params)):
            self.var.cond_params[i].smooth_factor = smooth_factor
            self.var.cond_params[i].k_cond = k_cond
            self.var.cond_params[i].k_evap = k_evap

    def set_particle_radii(self, radii):
        particle_radius = self.var.particle_radius
        for key in radii:
            ind = self.dat.species_names.index(key)
            particle_radius[ind,:] = radii[key]
        self.var.particle_radius = particle_radius
        self.update_vertical_grid(TOA_alt=self.var.top_atmos)

    def initialize_to_zT(self, z, T, Kzz, mix, P_surf):

        # Copy everything
        z, T, Kzz, mix = deepcopy(z), deepcopy(T), deepcopy(Kzz), deepcopy(mix)

        # Ensure mix sums to 1
        ftot = np.zeros(z.shape[0])
        for key in mix:
            ftot += mix[key]
        for key in mix:
            mix[key] = mix[key]/ftot

        # Compute mubar at all heights
        mu = {}
        for i,sp in enumerate(self.dat.species_names[:-2]):
            mu[sp] = self.dat.species_mass[i]
        mubar = np.zeros(z.shape[0])
        for key in mix:
            mubar += mix[key]*mu[key]

        # Compute pressure from hydrostatic equation
        P = compute_pressure_of_ZT(z, T, mubar, self.dat.planet_radius, self.dat.planet_mass, P_surf)

        # Calculate the photochemical grid
        z_top = z[-1]
        z_bottom = 0.0
        dz = (z_top - z_bottom)/self.var.nz
        z_p = np.empty(self.var.nz)
        z_p[0] = dz/2.0
        for i in range(1,self.var.nz):
            z_p[i] = z_p[i-1] + dz

        # Now, we interpolate all values to the photochemical grid
        P_p = 10.0**np.interp(z_p, z, np.log10(P))
        T_p = np.interp(z_p, z, T)
        Kzz_p = 10.0**np.interp(z_p, z, np.log10(Kzz))
        mix_p = {}
        for sp in mix:
            mix_p[sp] = 10.0**np.interp(z_p, z, np.log10(mix[sp]))
        k_boltz = const.k*1e7
        den_p = P_p/(k_boltz*T_p)

        # Update photochemical model grid
        self.update_vertical_grid(TOA_alt=z_top) # this will update gravity for new planet radius
        self.set_temperature(T_p)
        self.var.edd = Kzz_p
        usol = np.ones(self.wrk.usol.shape)*1e-40
        species_names = self.dat.species_names[:(-2-self.dat.nsl)]
        for sp in mix_p:
            if sp in species_names:
                ind = species_names.index(sp)
                usol[ind,:] = mix_p[sp]*den_p
        self.wrk.usol = usol

        # prep the atmosphere
        self.prep_atmosphere(self.wrk.usol)

    def initialize_robust_stepper(self, usol):
        rdat = self.rdat
        rdat.nsteps_total = 0
        rdat.nerrors_total = 0
        self.initialize_stepper(usol)
    
    def robust_step(self):

        rdat = self.rdat

        converged = False
        give_up = False

        if rdat.nsteps_total is None:
            raise PhotoException("You must first initialize a robust stepper with 'initialize_robust_stepper'")
        if self.var.nsteps_before_conv_check >= self.var.nsteps_before_reinit:
            raise PhotoException("`nsteps_before_conv_check` should be < `nsteps_before_reinit`")

        for i in range(1):
            try:
                tn = self.step()
                rdat.nsteps_total += 1
            except PhotoException as e:
                self.initialize_stepper(np.clip(self.wrk.usol.copy(),a_min=1.0e-40,a_max=np.inf))
                if rdat.nerrors_total > rdat.nerrors_before_giveup:
                    give_up = True
                    break

            # If converged, then return
            if self.wrk.tn > self.var.equilibrium_time:
                converged = True
                break

            # If converged, then return
            converged = self.check_for_convergence()
            if converged:
                break
            
            # Reinit if time to do that
            if self.wrk.nsteps > self.var.nsteps_before_reinit:
                self.initialize_stepper(np.clip(self.wrk.usol.copy(),a_min=1.0e-40,a_max=np.inf))
                break

            # Reinit if negative numbers
            if np.min(self.wrk.mix_history[:,:,0]) < rdat.min_mix_reset:
                self.initialize_stepper(np.clip(self.wrk.usol.copy(),a_min=1.0e-40,a_max=np.inf))
                break

        return give_up, converged

ATMOSPHERE_INIT = \
"""alt      den        temp       eddy                       
0.0      1          1000       1e6              
1.0e3    1          1000       1e6         
"""

@nb.njit()
def gravity(radius, mass, z):
    G_grav = const.G
    grav = G_grav * (mass/1.0e3) / ((radius + z)/1.0e2)**2.0
    grav = grav*1.0e2 # convert to cgs
    return grav

@nb.njit()
def hydrostatic_equation_z(z, u, planet_radius, planet_mass, ptm):
    P = u[0]
    grav = gravity(planet_radius, planet_mass, z)
    T, mubar = ptm.temperature_mubar(z)
    k_boltz = const.Boltzmann*1e7
    dP_dz = -(mubar*grav*P)/(k_boltz*T*const.Avogadro)
    return np.array([dP_dz])

@nb.experimental.jitclass()
class TempAltMubar:

    z : types.double[:] # type: ignore
    T : types.double[:] # type: ignore
    mubar : types.double[:] # type: ignore

    def __init__(self, z, T, mubar):
        self.z = z.copy()
        self.T = T.copy()
        self.mubar = mubar.copy()

    def temperature_mubar(self, z):
        T = np.interp(z, self.z, self.T)
        mubar = np.interp(z, self.z, self.mubar)
        return T, mubar

def compute_pressure_of_ZT(z, T, mubar, planet_radius, planet_mass, P_surf):

    ptm = TempAltMubar(z, T, mubar)
    args = (planet_radius, planet_mass, ptm)

    # Integrate to TOA
    out = integrate.solve_ivp(hydrostatic_equation_z, [z[0], z[-1]], np.array([P_surf]), t_eval=z, args=args, rtol=1e-6)
    assert out.success

    # Stitch together
    P = out.y[0]

    return P