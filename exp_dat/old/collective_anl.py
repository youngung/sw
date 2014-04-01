### In case that the 2th-th scan data files are provided collective in a single file...
import numpy as np
import os
from glob import glob
import matplotlib.pyplot as plt
from anl import inte

def read_th2(fn='02SEP13/eqa_mxk60_n7_10pct-1.uxd'):
    theta2, counts, dum1, dum2 = reader(fn)

    return theta2[0]

def reader(fn='02SEP13/eqa_mxk60_n7_10pct-1.uxd'):
    datstring = open(fn,'r').read()
    datb = datstring.split('; Data for range')
    header = datb[0]
    datb = datb[1:]
    print '# of samplings is ', len(datb)
    THETA2 = []
    COUNTS = []
    PHI = []
    KHI = []
    X, Y, Z = [], [], []

    for i in range(len(datb)):
        db = datb[i]
        dl = db.split('\n') # lines

        n = dl[0] # data number

        steptime = float(dl[2].split('=')[1])
        stepsize = float(dl[3].split('=')[1])

        th2_0 = float(dl[5].split('=')[1]) # included.
        th2_1 = float(dl[6].split('=')[1])

        khi = float(dl[8].split('=')[1])
        phi = float(dl[9].split('=')[1])
        KHI.append(khi); PHI.append(phi)

        x = float(dl[10].split('=')[1])
        y = float(dl[11].split('=')[1])
        z = float(dl[12].split('=')[1])

        X.append(x); Y.append(y); Z.append(z);

        dls = dl[50:]

        counts = []; theta2 = []
        for j in range(len(dls)-2):
            th2 = float(dls[j].split()[0])
            cnt = float(dls[j].split()[1]) / steptime
            counts.append(cnt)
            theta2.append(th2)

        COUNTS.append(counts)
        THETA2.append(theta2)

    COUNTS = np.array(COUNTS)
    THETA2 = np.array(THETA2)
    return THETA2,COUNTS,PHI,KHI

class Peak:
    def __init__(self,th2,cn,phi,khi,phase,hkl=None,TH2=None,dfc=False):
        self.th2=th2
        self.cn=cn
        self.phi=phi
        self.khi=khi
        self.dfc=dfc
        self.inte()
        self.phase = phase
        self.TH2 = np.array(self.th2).mean()
        self.hkl = hkl
    def inte(self):
        self.area = inte(self.th2, self.cn)
        if self.dfc:
            dfc=Defocus(fn='ass_0.dfc')
            self.area = self.area * dfc.calc_dfc(self.khi)


class Phase:
    def __init__(self,name):
        self.peaks = []
        self.name = name
        self.summed_peaks = []
        self.NP = 0

    def add_peak(self,peak):
        self.peaks.append(peak)

    def sumall(self,hkl=['','','']):
        """ Calculate average (summed) intensity profiles """
        from anl import bg_subtract as bgst
        TH2s = []
        theta2_axis = []
        for i in range(len(self.peaks)): # individual peaks
            dum = self.peaks[i].TH2      # its central 2theta angle
            if dum not in TH2s:
                TH2s.append(dum)
                theta2_axis.append(self.peaks[i].th2) # two thetas

        n2th = len(TH2s) # of unique 2theta
        self.NP = n2th

        print '%i peaks are probed'%n2th

        counts = []
        th2s = []
        for i in range(len(TH2s)): # counts for the 3 peaks in case of the austenite 2 peak in martensites
            counts.append([])
            th2s.append([])

        #th2s   = counts[::]
        tcounts = counts[::]

        for i in range(len(self.peaks)):
            TH2 = self.peaks[i].TH2 # central th2
            ind = np.searchsorted(TH2s,TH2) # to which unique peak, does this belong?
            dum = self.peaks[i].cn[::]

            counts[ind].append(dum)
            dum = self.peaks[i].th2[::]
            th2s[ind].append(dum)

        counts = np.array(counts)

        for i in range(len(counts)): # for 2 phases
            total_counts = np.zeros((len(counts[i][0]))) # for the individual scans
            for j in range(len(counts[i])): # individual scans

                nmulti = 1.0
                if self.peaks[j].khi==0:
                    nmulti = nmulti/2.
                if self.peaks[j].phi==0:
                    nmulti = nmulti/2.

                dum = counts[i][j] # individual scan # background subtraction

                bg_subtracted_cnts = bgst(theta2_axis[i], dum)

                for k in range(len(total_counts)): # along the 2theta
                    total_counts[k] = total_counts[k] + bg_subtracted_cnts[k] * nmulti

            total_counts = total_counts/len(counts[i][0])
            tcounts[i] = total_counts

        self.summed_theta2 = theta2_axis
        self.summed_counts = tcounts
        self.area=np.zeros((n2th))

        self.upeaks = []
        for i in range(n2th):
            dum = Peak(th2=self.summed_theta2[i], cn=self.summed_counts[i],
                       phi=-1,khi=-1,hkl=hkl[i],phase=self.name,dfc=True)
            self.upeaks.append(dum)

    def plot(self,color='k'):
        for ipeak in range(self.NP):
            plot(self.summed_theta[ipeak], self.summed_counts[ipeak])

class Scan:
    """
    fn='02SEP13/eqa_mxk60_n7_10pct-1.uxd',
    chosen_khiphi_ind=[],ieps_ph=True):
    """
    def __init__(self,fn='02SEP13/eqa_mxk60_n7_10pct-1.uxd',
                 chosen_khiphi_ind=[],ieps_ph=True):
        print self.__doc__

        self.ieps_ph = ieps_ph

        theta2, counts, PHI, KHI = reader(fn=fn)
        # gamma
        self.th2_axis = read_th2(fn)
        self.cnt = counts
        self.th2 = theta2
        self.fn = fn
        self.nscan = len(counts)
        print '# of scans:', self.nscan
        self.phis = PHI
        self.khis = KHI

        # self.ph
        self.ph = []
        # gamma
        self.gamma_hkl = ['200','220','311']
        self.gamma_2th = [59.67,89.63,111.38]
        # alpha
        self.alpha_hkl = ['200','211']#,'220']
        self.alpha_2th = [76.93,99.4]#,119.]]
        # epsilon?
        if ieps_ph: 
            self.epsilon_hkl = ['100']
            self.epsilon_2th   = [43.57]

        for i in range(len(self.gamma_hkl)):
            self.ph.append('gamma')
        for i in range(len(self.alpha_hkl)):
            self.ph.append('alpha')
        if ieps_ph: 
            for i in range(len(self.epsilon_hkl)):
                self.ph.append('epsilon')

        self.peaks_2th = self.gamma_2th + self.alpha_2th 
        if ieps_ph: self.peaks_2th = self.peaks_2th + self.epsilon_2th
        self.hkls  = self.gamma_hkl + self.alpha_hkl 
        if ieps_ph: self.hkls = self.hkls + self.epsilon_hkl
        self.npeaks = len(self.peaks_2th)


        # analsys conditions
        self.dt2 = 3.0

        # peaks
        self.chosen_khiphi = chosen_khiphi_ind
        self.peak_int()

        # phase analysis
        self.phase_assign()
        for iph in range(self.nph):
            if iph==0            : hkl = self.gamma_hkl
            if iph==1            : hkl = self.alpha_hkl
            if iph==2 and ieps_ph: hkl = self.epsilon_hkl
            self.phases[iph].sumall(hkl=hkl)

        # reflections
        self.reflections()
        print 'There are %i phases from this scan'%self.nph

    def reflections(self):
        reflcs = []
        for i in range(len(self.peaks_2th)):
            th2 = self.peaks_2th[i]
            th2 = round(th2,0)
            reflcs.append(th2)

        hkls   = self.hkls
        self.reflections = []
        for ir in range(len(reflcs)):
            self.reflections.append(
                Reflection(hkl=self.hkls[ir],th2=self.peaks_2th[ir],
                           phase=self.ph[ir]))

        for ip in range(len(self.peaks)):
            th2 = self.peaks[ip].TH2
            ind = reflcs.index(round(th2,0))
            self.reflections[ind].add_peak(self.peaks[ip])

    def peak_int(self):
        from anl import bg_subtract as bgst
        from anl import counts_2th as pick_a_peak

        self.peaks = []

        for i in range(self.nscan):
            if i in self.chosen_khiphi or len(self.chosen_khiphi)==0:
                c, t = self.cnt[i], self.th2[i]
                for j in range(self.npeaks):
                    tp, cp = pick_a_peak(self.peaks_2th[j], self.dt2, dat=[t,c])
                    cp = bgst(tp, cp)
                    apeak = Peak(th2=tp, cn=cp,
                                 phi=self.phis[i], khi=self.khis[i],
                                 phase=self.ph[j],dfc=False) # class initiation
                    self.peaks.append(apeak)

    def peak_sum(self):
        """ simple summation of all 2theta scans available. """
        for i in range(len(self.cnt)):
            if i==0: mn = 100000
            dum = len(self.cnt[i])
            if (dum<mn): mn = dum

        for i in range(len(self.cnt)):
            cnt = np.array(self.cnt[i])[:mn]
            if i==0: cnt_tot = cnt
            else: cnt_tot = cnt_tot + cnt
            th2 = self.th2[i][:mn]

            #plt.plot(cnt)

        self.cnt_tot = cnt_tot/len(self.cnt)

        plt.plot(th2, cnt_tot)

    def phase_assign(self):
        #ph = np.unique(self.ph)
        ph = ['gamma','alpha']#,'epsilon']
        if self.ieps_ph: ph.append('epsilon')
        self.nph = len(ph)
        self.phases = []
        for iph in range(self.nph):
            self.phases.append(Phase(name=ph[iph]))

        for i in range(len(self.peaks)): # 37(phi,khi) * 5 peaks = 185
            for iph in range(self.nph): #
                if self.peaks[i].phase==ph[iph]:
                    self.phases[iph].add_peak(self.peaks[i])

    def fnout(self):
        fn = self.fn.split('_')[-1].split('.')[0]
        fn = '%s.2th'%fn
        np.savetxt(fn,np.array([self.th2_axis, self.single_counts]).T)
        print '%s has been saved.'%fn

    def peak_analysis(self):
        from anl import Rhkl as RHKL

        self.gamma_n = len(self.gamma_hkl)
        self.alpha_n = len(self.alpha_hkl)
        if self.ieps_ph: self.epsilon_n = len(self.epsilon_hkl)

        gamma_area = []
        for ip in range(len(self.phases[0].upeaks)):
            dum = self.phases[0].upeaks[ip].area
            gamma_area.append(dum)

        alpha_area = []
        for ip in range(len(self.phases[1].upeaks)):
            dum = self.phases[1].upeaks[ip].area
            alpha_area.append(dum)

        if self.ieps_ph:
            epsilon_area = []
            for ip in range(len(self.phases[2].upeaks)):
                dum = self.phases[2].upeaks[ip].area
                epsilon_area.append(dum)

        self.gamma_area = gamma_area
        self.alpha_area = alpha_area
        if self.ieps_ph: self.epsilon_area = epsilon_area

        # R_HKL
        self.gamma_R = []
        for ip in range(len(self.phases[0].upeaks)):
            R = RHKL(hkl=self.gamma_hkl[ip],csym='fcc',th=self.gamma_2th[ip]/2.)
            self.gamma_R.append(R)

        self.alpha_R = []
        for ip in range(len(self.phases[1].upeaks)):
            R = RHKL(hkl=self.alpha_hkl[ip],csym='bcc',th=self.alpha_2th[ip]/2.)
            self.alpha_R.append(R)

        if self.ieps_ph:
            self.epsilon_R = []
            for ip in range(len(self.phases[2].upeaks)):
                R = RHKL(hkl=self.epsilon_hkl[ip],csym='hcp',th=self.epsilon_2th[ip]/2.)
                self.epsilon_R.append(R)

        ## volume fraction estimation
        gam_f = 0.
        for i in range(self.gamma_n):
            gam_f = gam_f + gamma_area[i]/ self.gamma_R[i]

        gam_f = gam_f / self.gamma_n

        alp_f = 0.
        for i in range(self.alpha_n):
            alp_f = alp_f + alpha_area[i]/ self.alpha_R[i]
        alp_f = alp_f / self.alpha_n

        if self.ieps_ph:
            eps_f = 0.
            for i in range(self.epsilon_n):
                eps_f = eps_f + epsilon_area[i]/ self.epsilon_R[i]

        f_sum = gam_f + alp_f
        if self.ieps_ph: f_sum = f_sum + eps_f
        v_gam = gam_f/f_sum
        v_alp = alp_f/f_sum
        if self.ieps_ph: v_eps = eps_f/f_sum

        self.vol = [v_gam, v_alp]#, v_eps]
        if self.ieps_ph: self.vol.append(v_eps)

class Reflection:
    def __init__(self,hkl,th2,phase):
        self.peaks = []
        self.hkl = hkl
        self.th2 = th2
        self.phase = phase
        pass
    def add_peak(self,peak):
        self.peaks.append(peak)

class Defocus:
    def __init__(self,fn='ass_0.dfc'):
        dat=np.loadtxt(fn,skiprows=5).T
        self.khi = dat[0]
        self.intensity = dat[1]
        self.dfc = dat[2]
    def calc_dfc(self,khi):
        return np.interp(khi,self.khi,self.dfc)

def vol_anal(paths='current',chosen_khiphi_ind=[],ieps_ph=True):
    import str_proj
    fns = glob('%s%s*.uxd'%(paths,os.sep))
    strains = []; Scans = []
    for i in range(len(fns)):
        eps = int(fns[i].split('pct')[0][-2:]) * 0.01
        strains.append(eps)
        myscan = Scan(fns[i],chosen_khiphi_ind=chosen_khiphi_ind,
                      ieps_ph=ieps_ph)
        myscan.peak_analysis()
        Scans.append(myscan)

    # mapping phi-khi sampling points using the equal area projection
    import eqa_proj
    fig3 = plt.figure(3)
    ax3  = fig3.add_subplot(121)
    ax4  = fig3.add_subplot(122)

    ax3.set_title('Equal area projection')
    ax4.set_title('Stereographic projection')
    eqa_proj.draw_circ(ax=ax3,nrim=0)
    eqa_proj.draw_circ(ax=ax4,nrim=0)

    Phi, Khi = [], []
    for i in range(len(Scans[0].peaks)):
        peak = Scans[0].peaks[i]
        p = peak.phi
        k = peak.khi
        Phi.append(p); Khi.append(k)
        # EQA projection
        x, y = eqa_proj.sph2cart(k*np.pi/180.,p*np.pi/180.)
        ax3.plot(x,y,'bx')
        # STR projection
        r, th = str_proj.sph2polar(p,k)   # unit: degree
        x, y  = eqa_proj.polar2cart(r,th*np.pi/180.) # unit: degree
        ax4.plot(x,y,'bx')

    ax3.set_ylim(-1.02,1.02);ax3.set_xlim(-1.02,1.02);ax3.set_axis_off()
    ax4.set_ylim(-1.02,1.02);ax4.set_xlim(-1.02,1.02);ax4.set_axis_off()

    grad = np.linspace(0.2,1.,len(strains))
    alpha = []; gamma = []; epsilon = []

    for i in range(len(Scans)):
        dat = Scans[i].vol
        if ieps_ph:
            gam, alp, eps = dat

            # if eps is negative?
            if eps<0: 
                gam = gam + eps
                eps = eps - eps # to fool the compiler

            if alp<0:
                gam = gam + alp
                alp = alp - alp # to fool the compiler

            if gam<0:
                alp = alp + gam
                gam = gam - gam # to fool the compiler            

        elif ieps_ph==False:
            gam, alp = dat

            if alp<0:
                gam = gam + alp
                alp = alp - alp # to fool the compiler            

            if gam<0:
                alp = alp + gam
                gam = gam - gam # to fool the compiler            

        alpha.append(alp)
        gamma.append(gam)
        if ieps_ph: epsilon.append(eps)

    fig1 = plt.figure(1); ax1 = plt.gca()
    ax1.plot(strains,gamma,'-+',label=r'$\gamma$')
    ax1.plot(strains,alpha,'-x',label=r'$\alpha$')
    if ieps_ph: ax1.plot(strains,epsilon,'-+',label=r'$\epsilon$')
    ax1.set_ylim(0.,1)
    ax1.set_ylabel(r'$f^{ph}$', dict(fontsize=28))
    ax1.set_xlabel(r'$\varepsilon$', dict(fontsize=28))

    fig2 = plt.figure(2); ax2 = fig2.add_subplot(111)
    colors = ['r','b','g']
    for i in range(len(Scans)):
        phs = Scans[i].phases
        for iph in range(len(phs)):
            for ip in range(len(phs[iph].upeaks)):
                peak = phs[iph].upeaks[ip]
                th2 = peak.th2
                cn = peak.cn
                #print th2, cn
                label = 'None'
                c = colors[iph]
                grd = grad[i]
                label =r'$\%s$'%phs[iph].name

                if ip==0 and i==len(Scans)-1:
                    ax2.plot(th2,cn,'-',color=c,label=label,alpha=grd)
                else: ax2.plot(th2,cn,'-',color=c,alpha=grd)

                x,y = th2[0],-150
                if i==0:
                    print phs[iph].upeaks[ip].hkl, phs[iph].name
                    tx = r'$(%s)_\%s$'%(phs[iph].upeaks[ip].hkl, phs[iph].name)
                    ax2.text(x,y,s=tx)
                     

    ## file for resulting volume fraction
    f = open('vf.txt', 'w')
    for i in range(len(Scans)):
        x  = strains[i]
        y0 = gamma[i]
        y1 = alpha[i]
        if ieps_ph: y2 = epsilon[i]
        f.writelines('%5.2f %7.4f %7.4f'%(x,y0,y1))
        if ieps_ph: f.writelines('%7.4f'%y2)
        f.writelines('\n')
    print 'Resulting VF %s has been saved.'%(f.name)
    f.close()

    ax2.set_xlabel(r'$2\theta$', dict(fontsize=28))
    ax2.set_ylabel('Intensity (counts)', dict(fontsize=28))
    ax1.legend(loc='best',fancybox=True).get_frame().set_alpha(0.5)
    ax2.legend(loc='best',fancybox=True).get_frame().set_alpha(0.5)

    plt.show()

    fig1.savefig('vol_basedon_EQA.pdf')
    fig2.savefig('peak_development.pdf')
    fig3.savefig('peaks_sampled.pdf')

    return Scans


def pfs(scans,ifig0=30,dfc=False,sym=False):
    for i in range(len(scans)):
        rot = 0
        if i==0: rot = 90.
        polefigure(scans[i],sym=sym,ifig=ifig0+i,dfc=True)

def polefigure(scan,sym=False,ifig=30,rot=0.,dfc=False):
    mydfc = Defocus('ass_0.dfc')
    import os
    #strain = float(scan.fn.split(os.sep)[1][:2])
    import str_proj, eqa_proj
    nph = scan.phases; hkls = scan.hkls; ths2 = []
    fig = plt.figure(ifig, figsize=(len(hkls)*1.0,0.65*len(hkls))); axes = []
    for i in range(len(scan.peaks_2th)):
        ths2.append(round(scan.peaks_2th[i],0))
        axes.append(fig.add_subplot(2,3,i+1))
        axes[i].set_axis_off()
        eqa_proj.draw_circ(ax=axes[i],nrim=0)
        axes[i].set_ylim(-1.1,1.1); axes[i].set_xlim(-1.1,1.1)
        axes[i].set_title(r'(%s)$_\%s$'%(
            scan.reflections[i].hkl,scan.reflections[i].phase))

    for i in range(len(scan.reflections)):
        mx = 0.; mn = 10.e10
        rfl = scan.reflections[i]
        peaks = rfl.peaks
        for j in range(len(peaks)):
            peak = peaks[j]
            area = peak.area
            if area>mx: mx = area
            if area<mn: mn = area
        for j in range(len(peaks)):
            k = peak.khi; p = peak.phi
            peak = peaks[j]
            area = peak.area
            # defocus correction
            if dfc: area = area * mydfc.calc_dfc(k)
            # STR projection
            r, th = str_proj.sph2polar(p,k)
            th = th + rot
            x, y  = eqa_proj.polar2cart(r,th*np.pi/180.)
            mn = mn * 0.90
            alpha = (area - mn) / (mx-mn) * 0.9
            ms = 6.
            axes[i].plot(x,y,'ko',ms=ms,mec='None',alpha=alpha)
            if sym:
                if x!=0:
                    x0 = -x
                    axes[i].plot(x0,y,'ko',ms=ms,mec='None',alpha=alpha)
                if y!=0:
                    y0 = -y
                    axes[i].plot(x,y0,'ko',ms=ms,mec='None',alpha=alpha)
                if y!=0 and y!=0:
                    x0 = -x
                    y0 = -y
                    axes[i].plot(x0,y0,'ko',ms=ms,mec='None',alpha=alpha)

    #print strain
    #fig.savefig('PF_%2ipct.pdf'%int(strain))

    # detect for each hkls.
