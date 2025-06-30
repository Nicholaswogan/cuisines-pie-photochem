import numpy as np
from photochem.clima import AdiabatClimate
import utils
import numba as nb

def test0_initialize(pc):
    "Initializes photochemical model to reasonable z-T-Kzz profile."

    c = AdiabatClimate(
        'slices/ArcheanEarth/test0/species_climate.yaml',
        'slices/ArcheanEarth/test0/settings_climate.yaml',
        'slices/ArcheanEarth/Sun_3.8Ga.txt'
    )
    P_i = np.ones(len(c.species_names))*1e-15
    P_i[c.species_names.index('H2O')] = 270
    P_i[c.species_names.index('N2')] = 0.7
    P_i[c.species_names.index('CO2')] = 0.3
    P_i *= 1e6
    c.T_trop = 180
    c.RH = np.ones(len(c.species_names))*0.5
    c.surface_temperature(P_i, 292)

    z = np.linspace(0,100e5,101)
    T = np.interp(z, c.z, c.T)
    P_surf = c.P_surf

    mix = {}
    species = ['H2O','N2','CO2']
    for i,sp in enumerate(species):
        ind = c.species_names.index(sp)
        mix[sp] = 10.0**np.interp(z, c.z, np.log10(c.f_i[:,ind]))

    z1, Kzz1 = np.loadtxt('slices/ArcheanEarth/test0/eddy_massie.txt',skiprows=2).T
    z1 = z1*1e5
    Kzz = 10.0**np.interp(z, z1, np.log10(Kzz1))

    pc.initialize_to_zT(z, T, Kzz, mix, P_surf)

    radii = {a: 1e-4 for a in pc.dat.species_names[:pc.dat.np]}
    radii['H2Oaer'] = 1e-3
    radii['HCaer1'] = 0.5e-4
    radii['HCaer2'] = 0.5e-4
    radii['HCaer3'] = 0.5e-4
    pc.set_particle_radii(radii)

def test0(savefile=False, use_atmosphere_file=True):

    atmosphere_file = None
    if use_atmosphere_file:
        atmosphere_file = 'slices/ArcheanEarth/test0/atmosphere.txt'

    pc = utils.EvoAtmosphereRobust(
        'slices/zahnle_earth_HNOCS.yaml',
        'slices/ArcheanEarth/test0/settings.yaml',
        'slices/ArcheanEarth/Sun_3.8Ga.txt',
        atmosphere_file
    )
    pc.set_particle_parameters(1, 1000, 10)

    if atmosphere_file is None:
        test0_initialize(pc)

    assert pc.find_steady_state()

    if savefile:
        pc.out2atmosphere_txt('slices/ArcheanEarth/test0/atmosphere.txt', overwrite=True)
        utils.pie_output_file(pc,'slices/ArcheanEarth/test0/Photochem_test0.txt')

def test1a_initialize(pc):
    z1, P1, T1, H2O1 = np.loadtxt('slices/ArcheanEarth/test1/ZPTH2O_Kasting1990_1km.txt',skiprows=1).T
    z1 *= 1e5
    z = np.linspace(0,100e5,201)
    T = np.interp(z, z1, T1)
    P_surf = P1[0]*1e6

    mix = {
        'N2': np.ones(len(z))*0.8,
        'CO2': np.ones(len(z))*0.2,
        'H2O': 10.0**np.interp(z, z1, np.log10(H2O1)),
    }

    z2, Kzz2 = np.loadtxt('slices/ArcheanEarth/test0/eddy_massie.txt',skiprows=2).T
    z2 *= 1e5
    Kzz = 10.0**np.interp(z, z2, np.log10(Kzz2))

    pc.initialize_to_zT(z, T, Kzz, mix, P_surf)

    radii = {a: 1e-4 for a in pc.dat.species_names[:pc.dat.np]}
    radii['H2Oaer'] = 1e-3
    radii['HCaer1'] = 0.5e-4
    radii['HCaer2'] = 0.5e-4
    radii['HCaer3'] = 0.5e-4
    pc.set_particle_radii(radii)

def test1a(savefile=False, use_atmosphere_file=True):

    atmosphere_file = None
    if use_atmosphere_file:
        atmosphere_file = 'slices/ArcheanEarth/test1/test1a/atmosphere.txt'

    pc = utils.EvoAtmosphereRobust(
        'slices/zahnle_earth_HNOC.yaml',
        'slices/ArcheanEarth/test1/test1a/settings.yaml',
        'slices/ArcheanEarth/Sun_3.8Ga.txt',
        atmosphere_file
    )
    pc.set_particle_parameters(1, 1000, 10)
    
    if atmosphere_file is None:
        test1a_initialize(pc)

    assert pc.find_steady_state()

    if savefile:
        pc.out2atmosphere_txt('slices/ArcheanEarth/test1/test1a/atmosphere.txt', overwrite=True)
        utils.pie_output_file(pc,'slices/ArcheanEarth/test1/test1a/Photochem_test1a.txt')

def test1d_180K(savefile=False, use_atmosphere_file=True):
    atmosphere_file = 'slices/ArcheanEarth/test1/test1a/atmosphere.txt'
    if use_atmosphere_file:
        atmosphere_file = 'slices/ArcheanEarth/test1/test1d/atmosphere_180K.txt'

    pc = utils.EvoAtmosphereRobust(
        'slices/zahnle_earth_HNOC.yaml',
        'slices/ArcheanEarth/test1/test1a/settings.yaml',
        'slices/ArcheanEarth/Sun_3.8Ga.txt',
        atmosphere_file
    )
    pc.set_particle_parameters(1, 1000, 10)

    # Temperature
    pc.set_temperature(np.ones(pc.var.nz)*180)

    assert pc.find_steady_state()

    if savefile:
        pc.out2atmosphere_txt('slices/ArcheanEarth/test1/test1d/atmosphere_180K.txt', overwrite=True)
        utils.pie_output_file(pc,'slices/ArcheanEarth/test1/test1d/Photochem_test1d_180K.txt')

def test1d_288K(savefile=False, use_atmosphere_file=True):
    # This one is hard to converge

    atmosphere_file = 'slices/ArcheanEarth/test1/test1d/atmosphere_288K.txt'

    pc = utils.EvoAtmosphereRobust(
        'slices/zahnle_earth_HNOC.yaml',
        'slices/ArcheanEarth/test1/test1a/settings.yaml',
        'slices/ArcheanEarth/Sun_3.8Ga.txt',
        atmosphere_file
    )
    pc.set_particle_parameters(1, 1000, 10)
    pc.var.atol = 1e-13

    assert pc.find_steady_state()

    if savefile:
        pc.out2atmosphere_txt('slices/ArcheanEarth/test1/test1d/atmosphere_288K.txt', overwrite=True)
        utils.pie_output_file(pc,'slices/ArcheanEarth/test1/test1d/Photochem_test1d_288K.txt')

def test2(savefile=False, use_atmosphere_file=True):

    atmosphere_file = None
    if use_atmosphere_file:
        atmosphere_file = 'slices/ArcheanEarth/test2/atmosphere.txt'

    pc = utils.EvoAtmosphereRobust(
        'slices/zahnle_earth_HNOC.yaml',
        'slices/ArcheanEarth/test2/settings.yaml',
        'slices/ArcheanEarth/test2/Sun_3.8Ga_energyunits_edited.txt',
        atmosphere_file
    )
    pc.set_particle_parameters(1, 1000, 10)
    
    if atmosphere_file is None:
        test1a_initialize(pc) # same as test1a

    assert pc.find_steady_state()

    if savefile:
        pc.out2atmosphere_txt('slices/ArcheanEarth/test2/atmosphere.txt', overwrite=True)
        utils.pie_output_file(pc,'slices/ArcheanEarth/test2/Photochem_test2.txt')

# def test3a(savefile=False):
#     print('\nRunning test3a')
#     pc = Atmosphere("reactions/zahnle_earth.yaml",\
#                 "./slices/ArcheanEarth/test3/test3a/settings_ArcheanEarth_test3a.yaml",\
#                 "./slices/ArcheanEarth/test2/Sun_3.8Ga_energyunits_edited.txt",\
#                 "./slices/ArcheanEarth/test3/test3a/atmosphere_ArcheanEarth_test3a.txt")
    
#     pc.initialize_stepper(pc.wrk.usol)
#     tn = 0.0
#     while tn < pc.var.equilibrium_time:
#         tn = pc.step()

#     if savefile:
#         utils.pie_output_file(pc,'slices/ArcheanEarth/test3/test3a/ArcheanEarth_test3a.txt')

# def test3b(savefile=False):
#     print('\nRunning test3b')
#     pc = Atmosphere("reactions/zahnle_earth.yaml",\
#                 "./slices/ArcheanEarth/test3/test3b/settings_ArcheanEarth_test3b.yaml",\
#                 "./slices/ArcheanEarth/test2/Sun_3.8Ga_energyunits_edited.txt",\
#                 "./slices/ArcheanEarth/test3/test3b/atmosphere_ArcheanEarth_test3b.txt")
    
#     pc.initialize_stepper(pc.wrk.usol)
#     tn = 0.0
#     while tn < pc.var.equilibrium_time:
#         tn = pc.step()

#     if savefile:
#         utils.pie_output_file(pc,'slices/ArcheanEarth/test3/test3b/ArcheanEarth_test3b.txt')

# def test3c(savefile=False):
#     print('\nRunning test3c')
#     pc = Atmosphere("reactions/zahnle_earth.yaml",\
#                 "./slices/ArcheanEarth/test3/test3c/settings_ArcheanEarth_test3c.yaml",\
#                 "./slices/ArcheanEarth/test2/Sun_3.8Ga_energyunits_edited.txt",\
#                 "./slices/ArcheanEarth/test3/test3c/atmosphere_ArcheanEarth_test3c.txt")
    
#     # Turn off molecular diffusion for all species but H2 and H
#     @nb.cfunc(nb.double(nb.double, nb.double, nb.double))
#     def custom_binary_diffusion_fcn(mu_i, mubar, T):
#         if mu_i < 2.01595:
#             b = 1.52e18*((1.0/mu_i+1.0/mubar)**0.5)*(T**0.5)
#         else:
#             b = 0.0    
#         return b
#     pc.var.custom_binary_diffusion_fcn = custom_binary_diffusion_fcn
    
#     pc.initialize_stepper(pc.wrk.usol)
#     tn = 0.0
#     while tn < pc.var.equilibrium_time:
#         tn = pc.step()

#     if savefile:
#         utils.pie_output_file(pc,'slices/ArcheanEarth/test3/test3c/ArcheanEarth_test3c.txt')

# def test3d(savefile=False):
#     print('\nRunning test3d')
#     pc = Atmosphere("reactions/zahnle_earth.yaml",\
#                 "./slices/ArcheanEarth/test3/test3d/settings_ArcheanEarth_test3d.yaml",\
#                 "./slices/ArcheanEarth/test2/Sun_3.8Ga_energyunits_edited.txt",\
#                 "./slices/ArcheanEarth/test3/test3d/atmosphere_ArcheanEarth_test3d.txt")
    
#     @nb.cfunc(nb.double(nb.double, nb.double, nb.double))
#     def custom_binary_diffusion_fcn(mu_i, mubar, T):
#         b = 1.52e18*((1.0/mu_i+1.0/mubar)**0.5)*(T**0.5)
#         return b
#     pc.var.custom_binary_diffusion_fcn = custom_binary_diffusion_fcn
    
#     pc.initialize_stepper(pc.wrk.usol)
#     tn = 0.0
#     while tn < pc.var.equilibrium_time:
#         tn = pc.step()

#     if savefile:
#         utils.pie_output_file(pc,'slices/ArcheanEarth/test3/test3d/ArcheanEarth_test3d.txt')

if __name__ == "__main__":
    savefile = True
    # test0(savefile=savefile, use_atmosphere_file=True)
    # test1a(savefile=savefile, use_atmosphere_file=True)
    # test1d_180K(savefile=savefile, use_atmosphere_file=False)
    # test1d_288K(savefile=savefile, use_atmosphere_file=False)
    test2(savefile=savefile, use_atmosphere_file=False)
    # test3a(savefile=savefile)
    # test3b(savefile=savefile)
    # test3c(savefile=savefile)
    # test3d(savefile=savefile)