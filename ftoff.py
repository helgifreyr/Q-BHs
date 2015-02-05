def funct_to_functf():
    functf = open('functf.dat','w')
    for line in open('funct.dat'):
        if len(line)>1:
            functf.write(line.split()[2]+'\n')

funct_to_functf()
