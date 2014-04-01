""" equal area projection """
import numpy as np
import matplotlib.pyplot as plt
cos = np.cos
sin = np.sin
sqrt = np.sqrt
acos = np.arccos
atan2 = np.arctan2
pi  = np.pi

def polar2sph(r,th):
    """
    polar coordinate to spherical coordinate (rot, tilting)

    arguemtns:
    r   (radius) * largest rmax=sqrt(2.)
    th  (planar rotation angle)
    """
    khi = 2. * np.arcsin(r/2.) # zenith (tilting khi)
    phi = th              # azimuth (planar rot)
    return khi, phi

def cart2polar(x,y):
    r = sqrt(x**2 + y**2)
    ph = atan2(y,x)
    return r, ph

def sph2polar(zenith, azimuth):
    """
    * zenith is the tilting 'khi' angle...
    sphearical cooridnates to polar coordinates (R, th)"""
    R = 2 * sin(zenith/2.) #/sqrt(2.)
    th = azimuth
    return R, th

def polar2cart(R,th):
    """ polar coordinate (R,th) to (x,y) """
    x = R * cos(th)
    y = R * sin(th)
    return x, y

def sph2cart(zenith,azimuth):
    """
    zenith (the tilting khi)
    azimuth (the plane rotation phi)
    """
    r, th = sph2polar(zenith=zenith, azimuth=azimuth) # r maximu is sqart(2.)
    r = r / sqrt(2.)
    x, y  = polar2cart(r, th)
    return x, y

def cart2sph(x,y):
    r, ph = cart2polar(x,y)
    r = r * sqrt(2.)
    khi, phi = polar2sph(r,ph) # khi: zenith, phi: azimuth
    return khi, phi

def check(max_khi = 90.,nkhi=9, nphi=15):
    max_khi = max_khi * pi/180.
    ax = plt.gca()
    draw_circ(ax)
    khi = np.linspace(0,max_khi,nkhi) # zenith
    phi = np.linspace(-pi,pi,nphi)  # azimuth

    print 'max khi', max(khi) * 180. / pi

    #khi = khi * pi / 180.
    #phi = phi * pi / 180.
    #khi1, phi1 = np.meshgrid(khi,phi)

    carts = []
    x0 = []
    for i in range(len(phi)):
        for j in range(len(khi)):
            p = phi[i]
            k = khi[j]

            x, y = sph2cart(zenith = k, azimuth=p)
            # kk, pp = cart2sph(x*sqrt(2.),y*sqrt(2.))
            # print 'k, kk', k*180./pi, kk*180./pi
            # print 'p, pp', p*180./pi, pp*180./pi
            # raw_input()
            carts.append([x,y])

            if i==0: x0.append(x)

    carts = np.array(carts)
    x, y = carts.T

    x = x / sqrt(2.)
    y = y / sqrt(2.)
    x0 = x0 / sqrt(2.)

    print 'max x', max(x)
    print 'x0:', x0
    print 'x0_diff', np.diff(x0)

    ax.plot(x,y,'.')
    ax.set_title('equal area projection')

def translate_diag(x,y,dx):
    x = x + dx/2.
    y = y + dx/2.*sqrt(3.)
    return x, y
    
def draw_circ(ax,nrim=4):
    phi = np.linspace(-pi,pi,1000)  # azimuth

    # multiplie rims
    chi = np.linspace(0,pi/2., nrim+2)[:-1:]
    for i in range(len(chi)):
        r, dum = sph2polar(chi[i], 0)
        r = r / sqrt(2.)
        xr, yr = r * cos(phi), r * sin(phi)    
        ax.plot(xr,yr,color='blue',ls='--',alpha=0.2)

    xr, yr = cos(phi), sin(phi)    
    ax.plot(xr,yr,color='b',ls='-')
    #ax.plot([0,0],'b.')
    ax.set_aspect('equal')

def main(khi_max=90., N=8):
    """ 
    """
    import str_proj as stereo
    khi_max = khi_max * pi / 180.
    fig1 = plt.figure(1); fig2 = plt.figure(2)
    ax1  = fig1.add_subplot(111); ax2 = fig2.add_subplot(111)
    draw_circ(ax1); draw_circ(ax2,nrim=0)
    x_mx, dum = sph2cart(khi_max,0.)
    y_mx = x_mx
    x = np.linspace(0.,x_mx,N)
    dx = np.diff(x)[0]
    dy = dx * sqrt(3)

    y_ori = 0.; y = []
    for i in range(len(x)):
        if y_ori>y_mx: break
        y.append(y_ori)
        y_ori = y_ori + dy
    y = np.array(y)
    x_pc, y_pc = [], []

    for i in range(len(x)):
        for j in range(len(y)):
            rp = sqrt(x[i]**2 + y[j]**2)
            if rp>x_mx: pass
            else:
                ax1.plot(x[i],y[j],'o',mfc='None', mec='r',
                         alpha=0.5)
                ax1.plot(x[i],y[j],'r.',alpha=0.5)
                x_pc.append(x[i]); y_pc.append(y[j])

            xt, yt = translate_diag(x[i],y[j],dx)
            rp = sqrt(xt**2+yt**2)
            if rp>x_mx: pass
            else:
                ax1.plot(xt,yt,'o',mfc='None', ms=0.1,mec='r',
                         alpha=0.5)
                ax1.plot(xt,yt,'r.',alpha=0.5)
                x_pc.append(xt); y_pc.append(yt)

    ## carts to polar
    khis, phis = [], []
    for i in range(len(x_pc)):
        x, y  = x_pc[i], y_pc[i]
        khi, phi = cart2sph(x,y)
        khis.append(khi); phis.append(phi)

    label = 'mxkhi%3s_N%s'%(
        str(int(round(khi_max*180./pi,0))).zfill(3),
        str(N).zfill(2))

    print 'label:',label

    f  = open('phi_khi_%s.txt'%label,'w')
    f1 = open('khi_phi_%s.txt'%label,'w')

    print '# of total angular points:', len(khis),\
        'Max khi:', khi_max*180.001/pi
    print 'It takes %3.1f hours provided that a '\
        'single scan takes 8 minutes'%(len(khis)*8./60.)
    f.writelines('%s\n'%label);  f1.writelines('%s\n'%label)
    print 'ind   phi   khi'
    f.writelines( 'ind   phi   khi\n');
    f1.writelines('ind   khi   phi\n')
    for i in range(len(khis)):
        khi, phi = khis[i], phis[i]
        print '%2i    %s  %s '%(
            i+1, str(round(phi*180./pi,1)).zfill(4),
            str(round(khi*180./pi,1)).zfill(4))
        f.writelines('%2i    %s  %s \n'%(
            i+1, str(round(phi*180./pi,1)).zfill(4),
            str(round(khi*180./pi,1)).zfill(4)))
        f1.writelines('%2i    %s  %s \n'%(
            i+1, str(round(khi*180./pi,1)).zfill(4),
            str(round(phi*180./pi,1)).zfill(4)))

        if i%5==4: 
            print '-----'
            f.writelines('-----\n')
            f1.writelines('-----\n')

        #raw_input()
        x, y = sph2cart(khi,phi)
        ax1.plot(x,y,'b+')
        # its corresponding stereographic projection
        sr, sth = stereo.sph2polar(
            azimuth=phi*180./pi, zenith=khi*180./pi)
        sth = sth * pi/180.
        sx, sy = polar2cart(sr,sth)
        #raw_input()
        ax2.plot(sx,sy,'o',mfc='None',mec='k',ms=3)

    print '\nphi chi angles are saved to %s'%f.name
    print '\nchi phi angles are saved to %s'%f1.name
    f.close(); f1.close()

    ax1.set_title(
        r'Equal Area projection $\chi^{max}:%i^\circ$'\
        '    N$^{tot}:%i$'%(
            int(khi_max*180.001/pi), N))
    ax2.set_title(r'Stereographic projection '\
                  '$\chi^{max}:%i^\circ$    N$^{tot}:%i$'%(
                      int(khi_max*180.001/pi), N))
    ax1.set_ylim(-1.1,1.1); ax1.set_xlim(-1.1,1.1)
    ax2.set_ylim(-1.1,1.1); ax2.set_xlim(-1.1,1.1)
    ax1.set_axis_off(); ax2.set_axis_off()
    fig1.savefig('eqa_proj_%s.pdf'%label)
    fig2.savefig('str_proj_%s.pdf'%label)
    #fig1.clf()

    
