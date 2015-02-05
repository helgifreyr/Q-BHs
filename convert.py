r_file = open('gridx.dat')
f_file = open('funct.dat')
new_f_file = open('new_funct.dat','w')

r_line = r_file.readline()
for f_line in f_file:
    if len(f_line)>10:
        new_f_line = f_line.split()
        new_f_line[-1] = str(float(f_line.split()[-1])*float(r_line))
        new_f_file.write(' '.join(new_f_line)+'\n')
    else:
        r_line = r_file.readline()
        new_f_file.write('\n')
