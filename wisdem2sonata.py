import yaml  
from jsonschema import validate
from SONATA.classAirfoil import Airfoil
from SONATA.classMaterial import read_IEA37_materials
from SONATA.classBlade import Blade
import numpy as np


with open('./jobs/PBortolotti/IEAonshoreWT.yaml', 'r') as myfile:
    inputs  = myfile.read()
with open('jobs/PBortolotti/IEAontology_schema.yaml', 'r') as myfile:
    schema  = myfile.read()
validate(yaml.load(inputs), yaml.load(schema))


yml = yaml.load(inputs)
airfoils = [Airfoil(af) for af in yml.get('airfoils')]
materials = read_IEA37_materials(yml.get('materials'))

byml = yml.get('components').get('blade')
B = Blade(name='TestBlade')
B.read_IEA37(byml, airfoils, materials, stations = [0.0, 0.1, 0.4], wt_flag = True)     

# B.iea37_converter(byml, airfoils, materials)

for key, cs in B.sections:
    print('STATUS:\t Building Section at grid location %s' % (key))
    cs.cbm_gen_topo()
    cs.cbm_gen_mesh()
    cs.cbm_run_vabs()
    title_plot = 'Blade station ' + str(key*100.) + '%'
    save_path  = 'jobs/PBortolotti/station_' + str(int(key*100.)) + '.pdf'
    cs.cbm_post_2dmesh(title=title_plot)