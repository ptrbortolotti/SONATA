

import subprocess


vabs_filename = 'anisopipe.dat'

cmd = ['wine ~/work/6_SONATA/SONATA/VABS/VABSIII.exe anisopipe.dat']

wine_dir = 'wine /Users/rfeil/work/6_SONATA/SONATA/VABS/VABSIII.exe'


vabs_in = ''

stdout = subprocess.run(cmd, stdout=subprocess.PIPE).stdout.decode('utf-8')

p = subprocess.call(['vim', abc])
