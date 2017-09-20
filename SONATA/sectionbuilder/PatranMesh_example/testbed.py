filename='SimpleMeshFile9.dat'
a = ''
with open(filename) as f:
    for line in f:
        line = line.partition('#')[0]
        line = line.rstrip()
        a += line 
        a += '\n'