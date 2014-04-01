import numpy as np
import matplotlib.pyplot as plt
from glob import glob
def reader(fn='0.30_full_2th_scan_phi00chi00.UXD'):
    dat = np.loadtxt(fn,skiprows=57).T
    th2 = dat[0]
    counts = dat[1]
    return th2, counts

def plotall(strain=0.3,ifig=2):
    fig=plt.figure(ifig); ax= fig.add_subplot(111)
    fns = glob('%3.2f_*.UXD'%strain)
    print fns
    for i in range(len(fns)):
        th2, counts = reader(fns[i])
        ax.plot(th2, counts)

def sumall(strain=0.3,ifig=1,path=None):
    import os
    if path==None: fn = '%3.2f_*.UXD'%strain
    else: fn = '%s%s%3.2f_*.UXD'%(path,os.sep,strain)
    fns = glob(fn)
    for i in range(len(fns)):
        th2, counts = reader(fns[i])
        if i==0: summed_counts = counts
        else: summed_counts = summed_counts + counts

    summed_counts = summed_counts / len(fns)

    if ifig!=None:
        fig=plt.figure(ifig); ax= fig.add_subplot(111)
        ax.plot(th2, summed_counts)

    fn = '%3.2f_total.2th'%strain
    np.savetxt(fn,np.array([th2,summed_counts]).T,fmt='%.3f   %i')
    print '%s has been saved.'%fn

def peak_anal(path=None):
    
    if path!=None:
        import os
        fns = glob('%s%s*.2th'%(path,os.sep))
    else:
        fns= glob('*.2th')
        
    fns.sort()
    strains = [0.]
    for i in range(len(fns)):
        strain = fns[i].split('_')[0]
        strains.append(strain)
    print fns
    raw_dat = []
    for i in range(len(fns)):
        print 'fns:', fns[i]
        raw_dat.append(np.loadtxt(fns[i]).T)

    d2th = 3. # for background subtraction
    ## Intensity based on area integration below
    gamma_hkl = ['200','220','311']
    gamma_2th = [59.67,89.63,111.38]

    alpha_hkl   = ['200','211']#,'220']
    alpha_2th = [76.93,99.4]#,119.]

    # for i in range(len(raw_dat)):
    #     deconvolution(raw_dat[i], 111, 51.5, ifig=i)
    # return

    vols = [[1,0.]]
    for i in range(len(raw_dat)):
        plt.plot(raw_dat[i][0], raw_dat[i][1],label=i)
        volf = vol(raw_dat[i], gamma_hkl, gamma_2th,
                   alpha_hkl, alpha_2th, d2th,'b')
        vols.append(volf)
        # volf = vol(raw_dat[1], gamma_hkl, gamma_2th,
        #            alpha_hkl, alpha_2th, d2th,'r')
        # vols.append(volf)

    vols = np.array(vols).T
    fig = plt.figure(33); ax = fig.add_subplot(111)
    fn = 'vf_exp.dat'
    fvol = open(fn, 'w')

    # strain phi1, phi2
    for i in range(len(strains)):
        fvol.writelines('%4.2f %.3f %.3f\n'%(
            float(strains[i]),vols[0][i], vols[1][i]))

    ax.plot(strains,vols[0],'o')
    ax.plot(strains,vols[1],'o')

    # fvol.writelines('%13s'%'ph\eps')
    # for i in range(len(strains)):
    #     fvol.writelines('%13s'%('%.3f'%float(strains[i])))
    # fvol.writelines('\n')
    # for i in range(len(vols)): # phases
    #     ax.plot(strains,vols[i],'o')
    #     fvol.writelines('%13s'%(' phase%i'%(i+1)))
    #     for j in range(len(vols[i])):
    #         fvol.writelines('%13s'%('%.3f'%vols[i][j]))
    #     fvol.writelines('\n')
    # fvol.close()

    print 'Volme fraction has been saved to %s.'%fn


    #ax.set_xlim(0.,); ax.set_ylim(0.,);


def vol(raw_dat,gamma_hkl, gamma_2th, alpha_hkl,
        alpha_2th, d2th,color='k'):
    from scipy.optimize import curve_fit
    gamma_n   = len(gamma_hkl)
    gamma_I   = []
    alpha_n   = len(alpha_hkl)
    alpha_I   = []
    for i in range(len(gamma_2th)):
        th2, cn = counts_2th(gamma_2th[i], d2th, raw_dat)
        cn = bg_subtract(th2, cn)
        gamma_I.append(inte(th2,cn))
        plt.plot(th2, cn, '-',color=color,
                 label=r'$(%s)_\gamma$'%gamma_hkl[i])

    for i in range(len(alpha_2th)):
        th2, cn = counts_2th(alpha_2th[i], d2th, raw_dat)
        cn = bg_subtract(th2, cn)
        INT_area = inte(th2,cn)
        # print 'Integrated peak:', INT_area
        alpha_I.append(INT_area)
        plt.plot(th2, cn, '--',color=color,
                 label=r'$(%s)_\alpha$'%alpha_hkl[i])

    ## Rhkl for each phases
    gamma_R = []
    #print 'gamma_R:'
    print '%5s %7s %11s'%('hkl','theta', 'R_hkl')
    for i in range(len(gamma_hkl)):
        R = Rhkl(hkl=gamma_hkl[i],csym='fcc',th = gamma_2th[i]/2.)
        print '%5s %7f %11e'%(gamma_hkl[i], gamma_2th[i]/2. ,R)
        gamma_R.append(R)

    alpha_R = []
    print '--'
    for i in range(len(alpha_hkl)):
        R = Rhkl(hkl=alpha_hkl[i],csym='bcc',th = alpha_2th[i]/2.)
        print '%5s %7f %11e'%(alpha_hkl[i], alpha_2th[i]/2. ,R)
        alpha_R.append(R)

    ## volume fraction estimation
    gam_f = 0.
    for i in range(gamma_n):
        gam_f = gam_f + gamma_I[i]/ gamma_R[i]
    gam_f = gam_f / gamma_n

    alp_f = 0.
    for i in range(alpha_n):
        alp_f = alp_f + alpha_I[i]/ alpha_R[i]
    alp_f = alp_f / alpha_n

    f_sum = gam_f + alp_f
    v_gam = gam_f/f_sum
    v_alp = alp_f/f_sum

    print 'vol %', v_gam*100., v_alp*100.

    plt.show()
    #plt.legend(loc='best')

    return v_gam, v_alp


def inte(th,cn):
    import scipy.integrate as integrate
    area = integrate.trapz(x=th,y=cn)
    return area

def bg_subtract(th,cn):
    x0 = th[0]
    x1 = th[-1]
    y0 = cn[0]
    y1 = cn[-1]

    slope = (y1-y0)/(x1-x0)

    counts = []
    for i in range(len(th)):
        x = th[i]
        y = __bg_func__(slope,x0,y0,x)
        newy = cn[i] - y
        #if newy<0: newy = 0.
        counts.append(newy)

    return counts


def __bg_func__(slope,x0,y0,x):
    return slope * (x-x0) + y0


def counts_2th(th2,d2th,dat):
    TH2, cnt = dat

    i0found = False
    i1found = False

    if TH2[-1]<th2+d2th or TH2[0]>th2-d2th:
        print TH2[0], TH2[-1]
        print th2-d2th, th2+d2th
        print 'd2th:', d2th
        raise IOError, 'reduce the d2th since the peak is not within.'


    for i in range(len(TH2)):
        if TH2[i] > th2-d2th and not(i0found):
            i0 = i
            i0found = True
        if TH2[i] > th2+d2th:
            i1 = i
            i1found = True
            break

    if not(i0found) or not(i1found): raise IOError, 'no both i0 and i1 found'

    #print 'i0, i1:', i0, i1

    theta2 = TH2[i0:i1:]
    counts = cnt[i0:i1:]
    return theta2, counts


def Rhkl(hkl,csym,th):
    """
    Rhkl for volume fraction estimation
    """
    # nano meter dimension
    a_a = 0.287 # alpha (bcc)
    a_g = 0.359 # gamma (fcc)
    a_e = 0.254 # epsilon (hcp) a axis
    c_e = 0.415 # epsilon (hcp) c axis
    csym0 = None

    a_a = a_a * 10**-9 # [m]
    a_g = a_g * 10**-9
    a_e = a_e * 10**-9
    c_e = c_e * 10**-9

    if csym=='bcc' or csym=='fcc':
        csym0 = 'cubic'
        if csym=='bcc': a = a_a
        if csym=='fcc': a = a_g

        # volume of the unit cell
        v = a**3

    elif csym=='hcp':
        csym0 = csym
        basal_area = 3./2. * np.sqrt(3) * a_e**2
        v = basal_area * c_e / 3.

    #print 'unit volume:', v

    # structure factor
    F = Fstrct(csym)
    #print 'structure factor: ', F
    # multiplicity factor
    p = phkl(hkl,csym0)
    #print 'Multiplicity factor:', p
    # The angular term
    ang_term = fbragg(th)
    #print 'The angular term:', ang_term
    # Temperature factor
    temp_term = temp(th)
    rst = (1./v**2) * (F**2) * p * ang_term * temp_term
    #print 'Rhkl:', rst
    return rst

def xtal_f(csym='bcc'):
    """ return crystal structure factor"""
    if csym=='bcc': return 2
    if csym=='fcc': return 3
    raise IOError, 'Unexpected crystal structure'

def Fstrct(csym='bcc',hkl='123'):
    if csym=='bcc': return 2
    if csym=='fcc': return 4
    if csym=='hcp': 
        h,k,l = hkl
    """
    0       when h+2k = 3n and l=odd
    1       when h+2k = 3n +- 1 and l= even
    sqrt(3) when h+2k = 3n +- 1 and l= odd
    sqrt(4) when h+2k = 3n and l=even
    """
    ieven = False
    if int(l)%2==1: ieven = False
    aux = int(h) + 2 * int(k)
    
    if aux%3==0 and not(ieven):
        return 0
    elif aux%3!=0 and ieven:
        return 1
    elif aux%3!=0 and not(ieven):
        return np.sqrt(3.)
    elif aux%3!=0 and ieven:
        return np.sqrt(4.)
    raise IOError, 'unexpected structure factor calculation'

def mu(a):
    """
    return unit cell volume fraction
    a_g = 0.359 # gamma
    a_a = 0.287 # alpha
    """
    return a**3

def fbragg(th):
    th = th * np.pi/180.
    rst = 1+np.cos(2*th)**2
    rst = rst / (np.sin(th)**2)/ np.cos(th)
    return rst


def temp(th):
    """
    """
    wave_length = 1.788 * 10e-10 # [m unit] for cobalt (bruker XRD in GIFT/POSTECH)
    rst = np.exp(- M(theta=th, wave_length=wave_length))
    return rst


def M(theta,wave_length):
    """ In case of Fe """
    wave_length = wave_length * (10**10) # wavelenght in angstrom
    debye_th = 430.  # Debye characteristic temperature of Fe
    temp = 296.
#    mass = 1.
    amu = 55.845
    avo = 6.022 * 10e+23
    h   = 6.626 * 10e-34
    k   = 1.381 * 10e-23 # Boltzmann's constant

    coeff = 1.15 * (10e4) * temp / amu / (debye_th**2)
#    ceff =

    x = debye_th  / temp
    phix = 0.694251

    # if x==1.4527: phix = 0.694251
    # else: raise IOError, 'find the proper value for x in Wolfram alpha'
    #print 'theta:', theta

    sinthlam = (np.sin(theta*np.pi/180.) / wave_length) ** 2

    #print '(sin(th)/lambda)^2 = ', sinthlam

    m  = coeff * (phix + x/4.) * sinthlam
    #print 'M:', m
    #print 'exp^(-M)', np.exp(-m)
    return m

def phkl(hkl='012',csym='cubic'):
    """
    Multiplicity factor for hkl
         FCC   p     BCC    p      HCP
         111   8     110    12
         200   6     200    6
         220   12    211    24
         311   24    220    12
         222   8     310    24
         400   6     222    8
         331  24     321    48
    # """
    # hkl = str(hkl)
    h, k, l = hkl
    h = int(h); k = int(k); l = int(l)


    nz = 0
    for i in range(3):
        #print 'hkl[i]', hkl[i]
        if int(hkl[i])==0:
            #print 'int(hkl[i])',int(hkl[i])
            nz = nz + 1

    nu = len(np.unique(hkl))

    if csym=='cubic':
        if nz==3: raise IOError, 'All zero hkl?'
        if nz==0:
            #print 'nu', nu, type(nu)
            if nu == 3: return 48
            if nu == 2: return 24
            if nu == 1: return 8
        if nz==1:
            nu = len(np.unique(hkl))
            if nu==3: return 24
            if nu==2: return 12
            raise IOError, 'Unexpected case found'
        if nz==2: return 6
        raise IOError, 'Unexpected end of the phkl'

    if csym=='hcp':
        hk = hkl[:2]
        nu_hk = len(np.unique(hk))
        nz_hk = 0
        for i in range(2):
            if int(hk[i])==0: nz_hk = nz_hk + 1

        if l==0:
            if nz_hk==1: return 6
            if nu_hk==2: return 12
            if nu_hk==1: return 6

        if l!=0:
            if nz_hk==2: return 2
            if nz_hk==1: return 12
            if nu_hk==1: return 12
            if nu_hk==2: return 24
        raise IOError, 'Unexpected end of the phkl'


# def phix(x):
#     """
#     use wolfram...

#     >>> int(y/e^y-1) dy from 0 to 1
#     """
#     print '(int(y/e^y-1) dy from 0 to x) * 1/x'
#     rst = raw_input(""" x is given as : """ )

#     return float(rst)



#  gaussina function fitting
def deconvolution(raw_dat, hkl, th2,ifig=None):
    from scipy.optimize import curve_fit

#    raw_dat = reader(fn='0.40_full_2th_scan_phi00chi15.UXD')
    # peak area, elution time, width of gaussian, exponentional damping
    # peak area, elution time, width of gaussian, exponentional damping
    parguess = [2000, 49.,  0.20, 0.5,#0.055,
                500, 52.5, 0.2, 0.5]#0.04]

    if ifig!=None: fig = plt.figure(ifig); ax=fig.add_subplot(111)

    th2, cn = counts_2th(th2, 5., raw_dat)
    cn = bg_subtract(th2, cn)# th2, cn

    ax.plot(th2, cn,'gray')
    fineth = np.linspace(th2[0], th2[-1],10000)
    ax.plot(fineth, two_peaks(fineth, *parguess),'b-', label='guess')


    popt, pcov = curve_fit(two_peaks, th2, cn, parguess)
    print 'popt:', popt

    ax.plot(th2, two_peaks(th2, *popt), 'rx')

    pars1 = popt[0:4]
    print 'pars1', pars1
    pars2 = popt[4:8]
    print 'pars2', pars2

    peak1 = asym_peak(fineth, pars1)
    peak2 = asym_peak(fineth, pars2)

    ax.plot(fineth, peak1,'g--', label='peak1')
    ax.plot(fineth, peak2,'g-', label='peak2')
    ax.legend(loc='best')
    

def asym_peak(t, pars):
    'from Anal. Chem. 1994, 66, 1294-1301'
    from scipy.special import erf
    a0 = pars[0]  # peak area
    a1 = pars[1]  # elution time
    a2 = pars[2]  # width of gaussian
    a3 = pars[3]  # exponential damping term
    f = (a0/2/a3*np.exp(a2**2/2.0/a3**2 + (a1 - t)/a3)
         *(erf((t-a1)/(np.sqrt(2.0)*a2) - a2/np.sqrt(2.0)/a3) + 1.0))
    return f

def sym_peak(t, pars):
    ## sigma = pars[0]
    ## mu = 
    f = 1/1

def two_peaks(t, *pars):    
    'function of two overlapping peaks'
    a10 = pars[0]  # peak area
    a11 = pars[1]  # elution time
    a12 = pars[2]  # width of gaussian
    a13 = pars[3]  # exponential damping term
    a20 = pars[4]  # peak area
    a21 = pars[5]  # elution time
    a22 = pars[6]  # width of gaussian
    a23 = pars[7]  # exponential damping term   
    p1 = asym_peak(t, [a10, a11, a12, a13])
    p2 = asym_peak(t, [a20, a21, a22, a23])
    return p1 + p2


