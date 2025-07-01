from photochem.utils import resave_mechanism_with_atoms
from photochem import zahnle_earth
from photochem.utils import stars
from photochem.utils._format import flowmap, blockseqtrue, yaml, MyDumper, Loader, mechanism_dict_with_atoms
from photochem.utils._convert_utils import generate_photo_yaml_entries, sort_photos

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

def FormatReactions_main(data):
    
    order = ['reverse-reactions','atoms','species','particles','reactions','missing']
    copy = data.copy()
    data.clear()
    for key in order:
        if key in copy.keys():
            data[key] = copy[key]
            
    # Atmos
    if 'atoms' in data:
        for i in range(len(data['atoms'])):
            data['atoms'][i] = flowmap(data['atoms'][i])
    
    # Species
    if 'species' in data:
        for i in range(len(data['species'])):
            
            if data['species'][i]['name'] == False:
                data['species'][i]['name'] = "NO"
            
            order = ['name', 'composition', 'condensate', 'thermo','note']
            copy = data['species'][i].copy()
            data['species'][i].clear()
            for key in order:
                if key in copy.keys():
                    data['species'][i][key] = copy[key]
                        
            data['species'][i]['composition'] = flowmap(data['species'][i]['composition'])
            if 'thermo' in data['species'][i].keys():
                
                order = ['model', 'reference-pressure','temperature-ranges','data']
                copy = data['species'][i]['thermo'].copy()
                data['species'][i]['thermo'].clear()
                for key in order:
                    if key in copy.keys():
                        data['species'][i]['thermo'][key] = copy[key]
                    
                data['species'][i]['thermo']['temperature-ranges'] = blockseqtrue(data['species'][i]['thermo']['temperature-ranges'])
                
                data['species'][i]['thermo']['data'] = [blockseqtrue(a) for a in blockseqtrue(data['species'][i]['thermo']['data'])]

    # Particles
    if 'particles' in data:
        for i in range(len(data['particles'])):
            data['particles'][i]['composition'] = flowmap(data['particles'][i]['composition'])
            if data['particles'][i]['formation'] == 'reaction':
                flowstyle = ['rate-constant','low-P-rate-constant','high-P-rate-constant','efficiencies']
                for key in flowstyle:
                    if key in data['particles'][i].keys():
                        data['particles'][i][key] = flowmap(data['particles'][i][key])
            elif data['particles'][i]['formation'] == 'saturation':
                flowstyle = ['parameters','vaporization','sublimation','super-critical']
                for key in flowstyle:
                    data['particles'][i]['saturation'][key] = flowmap(data['particles'][i]['saturation'][key])
            
    # Reactions
    if 'reactions' in data:
        for i in range(len(data['reactions'])):
            order = ['equation','type','rate-constant','rate-constants','low-P-rate-constant',
                     'high-P-rate-constant','duplicate','efficiencies','JPL','citation','rxn','note']
            copy = data['reactions'][i].copy()
            data['reactions'][i].clear()
            for key in order:
                if key in copy.keys():
                    data['reactions'][i][key] = copy[key]
                    
            flowstyle = ['rate-constant','low-P-rate-constant','high-P-rate-constant','efficiencies']
            for key in flowstyle:
                if key in data['reactions'][i].keys():
                    data['reactions'][i][key] = flowmap(data['reactions'][i][key])

            if 'rate-constants' in data['reactions'][i]:
                for j in range(len(data['reactions'][i]['rate-constants'])):
                    data['reactions'][i]['rate-constants'][j] = flowmap(data['reactions'][i]['rate-constants'][j])
                    
    return data

def test4_reactions():
    with open('slices/ArcheanEarth/test4/original_files/ArcheanHaze.yaml','r') as f:
        dat = yaml.load(f,Loader)

    def convert_sp_names(sp):
        if sp == 'CH23':
            sp = 'CH2'
        return sp

    for i in range(len(dat['species'])):
        dat['species'][i]['name'] = convert_sp_names(dat['species'][i]['name'])

    for i in range(len(dat['reactions'])):
        eq = dat['reactions'][i]['equation']
        react = [a.strip() for a in eq.split('=>')[0].split('+')]
        prod = [a.strip() for a in eq.split('=>')[1].split('+')]
        react = [convert_sp_names(a) for a in react]
        prod = [convert_sp_names(a) for a in prod]
        eq1 = ' + '.join(react) + ' => ' + ' + '.join(prod)
        # if eq != eq1:
        #     print(eq)
        #     print(eq1)
        #     print()
        dat['reactions'][i]['equation'] = eq1
        
    dat = mechanism_dict_with_atoms(dat, ['H','O','N','C'])

    reactions = []
    photolysis = []
    for rx in dat['reactions']:
        if rx['type'] != 'photolysis':
            reactions.append(rx)
        else:
            photolysis.append(rx)
    dat['reactions'] = reactions

    possible_photos = generate_photo_yaml_entries([sp['name'] for sp in dat['species']])
    missing, not_missing = sort_photos(photolysis, possible_photos)

    dat['reactions'] += possible_photos
    dat['missing'] = missing

    dat = FormatReactions_main(dat)
    with open('slices/ArcheanEarth/test4/ArcheanHaze_PhotochemPhotolysis.yaml','w') as f:
        yaml.dump(dat,f,Dumper=MyDumper,sort_keys=False,width=70)

if __name__ == '__main__':
    # reaction_mechanisms()
    # create_stellar_fluxes()
    test4_reactions()

