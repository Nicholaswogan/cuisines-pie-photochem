import numpy as np
from photochem import Atmosphere, zahnle_earth
from photochem.clima import AdiabatClimate
from utils import PhotochemClima

def test0(savefile=False, verbose=0, overwrite=False):
    print('\nRunning test0')
    pc = Atmosphere(zahnle_earth,\
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
        pc.out2atmosphere_txt('slices/ArcheanEarth/test0/atmosphere_ArcheanEarth_test0.txt',overwrite=overwrite)

def test1a(savefile=False, overwrite=False):
    print('\nRunning test1a')
    pc = Atmosphere(zahnle_earth,\
                "./slices/ArcheanEarth/test1/test1a/settings_ArcheanEarth_test1a.yaml",\
                "./slices/ArcheanEarth/Sun_4.0Ga.txt",\
                "./slices/ArcheanEarth/test1/test1a/atmosphere_ArcheanEarth_test1a.txt")
    
    pc.initialize_stepper(pc.wrk.usol)
    tn = 0.0
    while tn < pc.var.equilibrium_time:
        tn = pc.step()

    if savefile:
        pc.out2atmosphere_txt('slices/ArcheanEarth/test1/test1a/atmosphere_ArcheanEarth_test1a.txt',overwrite=overwrite)
        
if __name__ == "__main__":
    test0()
    test1a()