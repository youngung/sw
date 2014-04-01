### hexagonal grid
import numpy as np
import matplotlib.pyplot as plt
pi = np.pi

def main(khi_mx = 60.,N=8):
    khi_mx = khi_mx  * pi / 180.
    D_mx = 2 * np.sin(khi_mx/2.)

    print 'D_mx: ',D_mx

    y = np.zeros((N))
    for j in range(N):
        y_dum = (np.sqrt(3.) * D_mx) / (2. * N ) * j
        if abs(y_dum)>D_mx: 
            print 'y_dum exceeds the D_max'
            print y_dum, D_mx
        else: y[j] = y_dum

    M = np.arange(-N+1,N)
    x = np.zeros((len(M),N))
    x[::] = -10
    khi = np.zeros((len(M),N))
    phi = np.zeros((len(M),N))
    D = np.zeros((len(M),N))

    for j in range(N):
        for k in range(len(M)):
            i = M[k]
            if j%2==0:
                dum = D_mx/N * i
            if j%2==1:
                dum = D_mx/(2*N) + D_mx/N * i

            if abs(dum)>np.sqrt(D_mx**2 - y[j]**2):
            #if abs(dum)>np.sqrt((2*np.cos(khi_mx/2.))**2 - y[j]**2):
                #raw_input()
                #print "exceeding!!"
                pass
            else:
                x[k,j] = dum
                D[k,j] = np.sqrt(x[k,j]**2 + y[j]**2)
                khi[k,j] = 2 * np.arcsin(D[k,j]/2.)
                if x[k,j]>0: 
                    phi[k,j] = np.arctan(y[j] / x[k,j])
                    #print phi[k,j]
                elif x[k,j]==0: phi[k,j] = pi / 2.
                elif x[k,j]<0: phi[k,j] = np.arctan(y[j]/x[k,j]) + pi

    # for i in range(len(x)):
    #     n = x[i]
    #     for j in range(len(n)):
    #         if x[i][j]==-10:
    #             pass
    #         else:
    #             if np.sqrt(x[i,j]**2+y[j]**2)>D_mx: pass
    #             else: plt.plot(x[i,j],y[j],'o', mfc='None')

    # plt.gca().set_ylim(0,1)
    # plt.gca().set_xlim(-1,1)

    for ix in range(len(phi)):
        for iy in range(len(phi[ix])):
            r, th = lambert_eqarea(khi[ix,iy], phi[ix,iy])
            x = r * np.cos(th)
            y = r * np.sin(th)

            plt.plot(x,y,'o',mfc='None')

    return phi, khi

def lambert_eqarea(khi,phi):
    """ equal area projection to spherical coorindate"""
    r = 2 * np.sin(khi/2.)
    th = phi
    return r, th
    
    
