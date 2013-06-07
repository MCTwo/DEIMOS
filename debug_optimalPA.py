import numpy
def objectPA(H,delta,phi):
    from math import atan, tan, cos, sin, pi
    d2r = pi/180.
    phi *= d2r    
    H = H*15*d2r
    delta *= d2r
    
    denom = tan(phi)*cos(delta)-sin(delta)*cos(H)
    q = atan(sin(H)/denom)
    q = q/d2r
    if denom < 0:
        q += 180
    return q

def optimalPA(H,delta,phi):
    '''
    This is based on:
    Filippenko, A.V., 1982. The importance of atmospheric differential refraction in spectrophotometry. Publications of the Astronomical Society of the Pacific, 94, pp.715–721. Available at: http://adsabs.harvard.edu/abs/1982PASP...94..715F.
    
    Input:
    phi = [float; units=degrees] observers latitude
    H = [float; units=hours] object's hour angle (H is + if west of the meridian)
    delta = [float; units=degrees] object's declination
    '''
    from math import pi,sin,cos,asin
    pa_obj = objectPA(H,delta,phi)
    d2r = pi/180.
    phi *= d2r
    if H < 0:
        sign = -1
        H = -H
    else:
        sign = 1
    H = H*15*d2r
    delta *= d2r
    eta_rad = sign*asin(sin(H)*cos(phi) / 
                        (1-(sin(phi)*sin(delta) +
                            cos(phi)*cos(delta)*cos(H))**2)**(0.5))
    eta_deg = eta_rad/d2r
    
    if sign*pa_obj < 0:
        eta_deg = sign*(180-sign*eta_deg)
    return eta_deg


dec = numpy.arange(-35,90,5)
HA = numpy.arange(0,8.5,0.5)
N_dec = numpy.size(dec)
N_HA = numpy.size(HA)
result = numpy.zeros((N_dec,N_HA))
pa = numpy.zeros((N_dec,N_HA))
for i in numpy.arange(N_dec):
    for j in numpy.arange(N_HA):
        result[i,j]=optimalPA(HA[j],dec[i],33.35)
        pa[i,j] = objectPA(HA[j],dec[i],33.35)
numpy.savetxt('palomar.txt',result,fmt='%1.1f',delimiter='\t')
numpy.savetxt('palomar_pa.txt',pa,fmt='%1.1f',delimiter='\t')
print result