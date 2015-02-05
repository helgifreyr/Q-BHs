from os import listdir
import glob
import scipy

# data format:
# 0   1     2   3   4   5     6    7   8    9     10     11   12   13   14  15
# rh, alfa, c1, c2, c3, Jint, THc, Mc, AHc, err1, minF0, f0H, f1H, f2H, ZH, Mint
# 16  17  18 19         20     21     22     23     24     25    26    27    28    29
# Le, Lp, w, constJINF, maxF0, minF1, maxF1, minF2, maxF2, minW, maxW, minZ, maxZ, errSmarr
def get_data(ms,alphas,ws,branch,rhs):
    string = 'm='+ms+'/alpha='+alphas+'/w='+ws+'/'+branch+'/rh='+rhs+'/'
    dirs = sorted(glob.glob(string))
    datas = scipy.zeros((len(dirs),30))
    i=0
    for dir in dirs:
        rep = ('Min',''), ('Max',''), (r'[',''), (r'],',''), (r'$Failed','-1'), (r'*^','e'), ('\n',''), (r'{',''), (r'}',''), ('EndOfFile','')
        data = ''.join([reduce(lambda a, kv: a.replace(*kv), rep, line) for line in open(dir+'/tmp.txt')]).split(',')
        print dir
        datas[i] = data
        i+=1
    return datas

def get_Bstar_data(ms,alphas,ws,branch):
    string = 'solitonic-limit/m='+ms+'/alpha='+alphas+'/'+branch+'/w='+ws+'/'
    dirs = sorted(glob.glob(string))
    datas = scipy.zeros((len(dirs),15))
    i=0
    for dir in dirs:
        rep = ('Min',''), ('Max',''), (r'[',''), (r'],',''), (r'$Failed','-1'), (r'*^','e'), ('\n',''), (r'{',''), (r'}',''), ('EndOfFile','')
        data = ''.join([reduce(lambda a, kv: a.replace(*kv), rep, line) for line in open(dir+'/tmp-py.txt')]).split()
        datas[i] = data
        i+=1
    return datas

def get_extremal_data(ms,alphas,ws,branch,rhs):
    string = 'extremal/m='+ms+'/alpha='+alphas+'/w='+ws+'/'+branch+'/rh='+rhs+'/'
    dirs = sorted(glob.glob(string))
    datas = scipy.zeros((len(dirs),30))
    i=0
    for dir in dirs:
        rep = ('Min',''), ('Max',''), (r'[',''), (r'],',''), (r'$Failed','-1'), (r'*^','e'), ('\n',''), (r'{',''), (r'}',''), ('EndOfFile','')
        data = ''.join([reduce(lambda a, kv: a.replace(*kv), rep, line) for line in open(dir+'/tmp.txt')]).split(',')
        datas[i] = data
        i+=1
    return datas
