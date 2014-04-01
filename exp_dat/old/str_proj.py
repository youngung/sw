## stereo projection
import numpy as np
sin = np.sin
cos = np.cos
pi = np.pi

def sph2polar(azimuth, zenith):
    """
    azimuth: rotation      (phi)
    zenith : tilting angle (khi)

    returns r, phi

    khi range: 180 ~ 90

    """
    khi = zenith * pi / 180. - pi
    phi = azimuth * pi / 180.

    if cos(khi)==1: r = 0
    else: r = sin(khi) / (cos(khi) - 1)

    # if sin(khi)==1: r = 0
    # else: r = cos(khi) / (1 - sin(khi))

    th = phi * 180. / pi
    #print 'r:', r, 'theta:', th

    return r, th

# def sph2xy(azimuth, zenith):
#     """
#     azimuth: rotation      (phi)
#     zenith : tilting angle (khi)

#     returns x,y

#     """
#     khi = zenith * pi / 180.
#     phi = azimuth * pi / 180.

#     x,y,z = spherical2cartesian(phi,khi)

#     X = x/(1-z)
#     Y = y/(1-z)

#     # if cos(khi) == 1: r = 0
#     # else: r = sin(khi) / (1 - cos(khi))

#     # # at khi = 90: r = 1
#     # # at khi = 0 : sin(khi) = 0, cos(khi) = 1: r = 0
#     # th = phi

#     return X,Y


# def spherical2cartesian(azimuth, zenith):
#     """
#     azimuth: rotation      (phi)
#     zenith : tilting angle (khi)
#     """
#     khi = zenith * pi / 180.
#     phi = azimuth * pi / 180.

#     x = sin(khi) * cos(phi)
#     y = sin(khi) * sin(phi)
#     z = cos(khi)

#     return x,y,z


    
