"""
"""
import numpy as np
from glob import glob
fns = glob('*.csv')
def reader(fn='fe110.csv'):
    dat=np.loadtxt(fn,skiprows=1,delimiter=',').T
    header= ['dis(mm)', 'load(kN)','Stress(MPa)',\
        'Xc','Xc_err','2theta','epshkl',\
        'epshkl_e','Width','Width_e','Area','Area_a','Height']

    # sigma11, epshkl, epshkl_e, Area, Area_a

    fnout = '%s.txt'%fn.split('.csv')[0]
    d = np.array([dat[2],dat[6],dat[7],dat[10],dat[11]])
    np.savetxt(fnout,d.T,fmt='%+11.4e')

def main():
    for i in range(len(fns)):
        reader(fns[i])
