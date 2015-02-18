pro RAdecfromPAdeltay, PA, RA, Dec, deltay

	if PA gt 0 then begin
		deltaRApix = (deltay/1.5126)*(sqrt(1+1/(tan(PA*2*!pi/360)^2)))^(-1)
		deltaDecpix = deltaRApix/tan(PA*2*!pi/360)
	end
	
	if PA lt 0 then begin
		deltaRApix = (deltay/1.5126)*sqrt((1+tan((90+PA)*2*!pi/360)^2))
		deltaDecpix = -deltaRApix*tan((90+PA)*2*!pi/360)
	end

	deltaDec = deltaDecpix * 0.184 * 1/3600
	deltaRA = deltaRApix * 0.184 * 1/3600/((cos(Dec*2*!pi/360)+cos((Dec+deltaDec)*2*!pi/360))/2)
	
	newDec = Dec + deltaDec
	newRA = RA + deltaRA
	
	print, newRA, newDec, format ='(f12.7, f11.7)'
end	
