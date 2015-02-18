pro PAdeltayfromRARAdecdec, RA, Dec, RA2, Dec2

	RArad = RA*2.*!PI/360.
        decrad = Dec*2.*!PI/360.
	RArad2 = RA2*2.*!PI/360.
        decrad2 = Dec2*2.*!PI/360.


	PA = 360/(2*!pi)*atan( (-cos(decrad)*sin(RArad-RArad2))/(cos(decrad)*sin(decrad2)-sin(decrad)*cos(decrad2)*cos(RArad - RArad2)))
       if PA lt 0 and RArad gt RArad2 then PA = 360. + PA else if PA lt 0 and RArad lt RArad2 then PA = 180. + PA else if PA gt 0 and RArad gt RArad2 then PA = 180. + PA else if PA gt 0 and RArad lt RArad then PA = PA
        d = 360/(2*!pi)*3600*acos(sin(decrad)*sin(decrad2)+cos(decrad)*cos(decrad2)*cos(RArad-RArad2))

	print, 'The separation between the two galaxies is', strcompress(string(d), /REMOVE_ALL), 'at a position angle of', strcompress(string(PA), /REMOVE_ALL)

end
