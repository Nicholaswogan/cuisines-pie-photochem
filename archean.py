import numpy as np
from photochem import Atmosphere
from photochem.clima import AdiabatClimate
from utils import PhotochemClima
import utils
import numba as nb

def test0(savefile=False, verbose=0):
    print('\nRunning test0')
    pc = Atmosphere("reactions/zahnle_earth.yaml",\
                "./slices/ArcheanEarth/test0/settings_ArcheanEarth_test0.yaml",\
                "./slices/ArcheanEarth/Sun_4.0Ga.txt",\
                "./slices/ArcheanEarth/test0/atmosphere_ArcheanEarth_test0.txt")
    pc.var.verbose=verbose
    
    c = AdiabatClimate('climate-settings/adiabat_species.yaml',\
                       'climate-settings/adiabat_settings.yaml',\
                       'slices/ArcheanEarth/Sun_4.0Ga.txt')
    c.T_trop = 200
    c.RH = np.ones(c.RH.shape[0])*1.0

    p = PhotochemClima(pc, c)
    p.iterate()

    if savefile:
        utils.pie_output_file(pc,'slices/ArcheanEarth/test0/ArcheanEarth_test0.txt')

def test1a(savefile=False):
    print('\nRunning test1a')
    pc = Atmosphere("reactions/zahnle_earth.yaml",\
                "./slices/ArcheanEarth/test1/test1a/settings_ArcheanEarth_test1a.yaml",\
                "./slices/ArcheanEarth/Sun_4.0Ga.txt",\
                "./slices/ArcheanEarth/test1/test1a/atmosphere_ArcheanEarth_test1a.txt")
    
    pc.initialize_stepper(pc.wrk.usol)
    tn = 0.0
    while tn < pc.var.equilibrium_time:
        tn = pc.step()

    if savefile:
        utils.pie_output_file(pc,'slices/ArcheanEarth/test1/test1a/ArcheanEarth_test1a.txt')

def test1d_180K(savefile=False):
    print('\nRunning test1d_180K')
    pc = Atmosphere("reactions/zahnle_earth.yaml",\
                "./slices/ArcheanEarth/test1/test1d/settings_ArcheanEarth_test1d_180K.yaml",\
                "./slices/ArcheanEarth/Sun_4.0Ga.txt",\
                "./slices/ArcheanEarth/test1/test1d/atmosphere_ArcheanEarth_test1d_180K.txt")
    
    pc.initialize_stepper(pc.wrk.usol)
    tn = 0.0
    while tn < pc.var.equilibrium_time:
        tn = pc.step()

    if savefile:
        utils.pie_output_file(pc,'slices/ArcheanEarth/test1/test1d/ArcheanEarth_test1d_180K.txt')

def test1d_288K(savefile=False):
    print('\nRunning test1d_288K')
    pc = Atmosphere("reactions/zahnle_earth.yaml",\
                "./slices/ArcheanEarth/test1/test1d/settings_ArcheanEarth_test1d_288K.yaml",\
                "./slices/ArcheanEarth/Sun_4.0Ga.txt",\
                "./slices/ArcheanEarth/test1/test1d/atmosphere_ArcheanEarth_test1d_288K.txt")
    
    pc.initialize_stepper(pc.wrk.usol)
    tn = 0.0
    while tn < pc.var.equilibrium_time:
        tn = pc.step()

    if savefile:
        utils.pie_output_file(pc,'slices/ArcheanEarth/test1/test1d/ArcheanEarth_test1d_288K.txt')

def test2(savefile=False):
    print('\nRunning test2')
    pc = Atmosphere("reactions/zahnle_earth.yaml",\
                "./slices/ArcheanEarth/test2/settings_ArcheanEarth_test2.yaml",\
                "./slices/ArcheanEarth/test2/Sun_3.8Ga_energyunits_edited.txt",\
                "./slices/ArcheanEarth/test2/atmosphere_ArcheanEarth_test2.txt")
    
    pc.initialize_stepper(pc.wrk.usol)
    tn = 0.0
    while tn < pc.var.equilibrium_time:
        tn = pc.step()

    if savefile:
        utils.pie_output_file(pc,'slices/ArcheanEarth/test2/ArcheanEarth_test2.txt')

def test3a(savefile=False):
    print('\nRunning test3a')
    pc = Atmosphere("reactions/zahnle_earth.yaml",\
                "./slices/ArcheanEarth/test3/test3a/settings_ArcheanEarth_test3a.yaml",\
                "./slices/ArcheanEarth/test2/Sun_3.8Ga_energyunits_edited.txt",\
                "./slices/ArcheanEarth/test3/test3a/atmosphere_ArcheanEarth_test3a.txt")
    
    pc.initialize_stepper(pc.wrk.usol)
    tn = 0.0
    while tn < pc.var.equilibrium_time:
        tn = pc.step()

    if savefile:
        utils.pie_output_file(pc,'slices/ArcheanEarth/test3/test3a/ArcheanEarth_test3a.txt')

def test3b(savefile=False):
    print('\nRunning test3b')
    pc = Atmosphere("reactions/zahnle_earth.yaml",\
                "./slices/ArcheanEarth/test3/test3b/settings_ArcheanEarth_test3b.yaml",\
                "./slices/ArcheanEarth/test2/Sun_3.8Ga_energyunits_edited.txt",\
                "./slices/ArcheanEarth/test3/test3b/atmosphere_ArcheanEarth_test3b.txt")
    
    pc.initialize_stepper(pc.wrk.usol)
    tn = 0.0
    while tn < pc.var.equilibrium_time:
        tn = pc.step()

    if savefile:
        utils.pie_output_file(pc,'slices/ArcheanEarth/test3/test3b/ArcheanEarth_test3b.txt')
        
if __name__ == "__main__":
    savefile = False
    test0(savefile=savefile)
    test1a(savefile=savefile)
    test1d_180K(savefile=savefile)
    test1d_288K(savefile=savefile)
    test2(savefile=savefile)
    test3a(savefile=savefile)
    test3b(savefile=savefile)