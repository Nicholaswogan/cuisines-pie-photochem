from photochem.utils import resave_mechanism_with_atoms
from photochem import zahnle_earth
from photochem.utils import stars

def reaction_mechanisms():
    resave_mechanism_with_atoms(
        zahnle_earth,
        'slices/zahnle_earth_HNOCS.yaml',
        ['H','O','N','C','S']
    )
    resave_mechanism_with_atoms(
        zahnle_earth,
        'slices/zahnle_earth_HNOC.yaml',
        ['H','O','N','C']
    )

def create_stellar_fluxes():
    _ = stars.solar_spectrum(
        outputfile='slices/ArcheanEarth/Sun_3.8Ga.txt',
        age=3.8,
        stellar_flux=1367,
    )

if __name__ == '__main__':
    reaction_mechanisms()
    create_stellar_fluxes()

