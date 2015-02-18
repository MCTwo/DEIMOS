; program to convert absolute magnitude to luminosity in units of ergs/s
function Ez, z

   distance = 1/sqrt(0.3*(1+z)^3+0.7)
   return, distance

end

pro MtoL, h, z1, M, restlam

Dh = 3000./h

mpctocm = double(3.08568e24)

Dcnorm = qromb('Ez',0.,z1)
DLnorm = Dcnorm*(1+z1)

Dc = Dcnorm*Dh
DL = DLnorm*Dh

Dccm = Dc*mpctocm
DLcm = DL*mpctocm

mapp = M + 5*alog10(DL*10.^6/10.)


fnu = 10.^(-(mapp+48.6)/2.5)

flam = fnu*3.e18/(restlam*(1.+z1))^2

Llam = 4.*!PI*DLcm^2*flam
L = Llam*restlam

print, 'The luminosity of a galaxy with absolute magnitude ', strcompress(string(M),/remove_all), ' observed at ', strcompress(string(restlam),/remove_all), ' Angstroms in the rest-frame is ', strcompress(string(L),/remove_all), ' ergs/s'

end 
