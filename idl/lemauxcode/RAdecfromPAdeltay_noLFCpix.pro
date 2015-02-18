pro RAdecfromPAdeltay_noLFCpix, PA, RA, Dec, deltay

	r = deltay
	if PA gt 0 then begin 
		deltadec = 0.1185*r*cos(PA*2*!pi/360)
		deltaRA = 0.1185*r*sin(PA*2*!pi/360)/sqrt(cos(Dec*2*!pi/360)*cos((Dec+deltadec)*2*!pi/360))
	end

	if PA lt 0 then begin 
		deltadec = -0.1185*r*sin((90+PA)*2*!pi/360)
		deltaRA = 0.1185*r*cos((90+PA)*2*!pi/360)/sqrt(cos(Dec*2*!pi/360)*cos((Dec+deltadec)*2*!pi/360))
	end

	newDec = Dec + deltaDec/3600
	newRA = RA + deltaRA/3600
	
	a = 99
        print, newRA, newDec, a, a, a, format ='(f12.7, 4x, f11.7, 6x, f7.4, 1x,  f7.4, 1x, f7.4)'
end
		
