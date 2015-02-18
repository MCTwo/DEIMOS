pro finderfromLFC5281, setup, maketable, viewchart, makechart

	hdr = headfits('/Volumes/Data2/orelse/lemaux/LFC/N5281/nep5281_deep_norm.fits')

	CRVAL1 = double(sxpar(hdr,'CRVAL1'))	;initial RA for any movement along the image
	CRVAL2 = double(sxpar(hdr, 'CRVAL2'))	;initial dec for any movement along the image
	
	aa = double(sxpar(hdr, 'CD1_1'))	;get the transformation matrix between pixels and global coordinates 
	ab = double(sxpar(hdr, 'CD1_2'))
	ba = double(sxpar(hdr, 'CD2_1'))
	bb = double(sxpar(hdr, 'CD2_2'))

	refxpix = sxpar(hdr, 'CRPIX1')		;initial x pix corresponding to CRVAL1
	refypix = sxpar(hdr, 'CRPIX2')		;initial y pix coresponding to CRVAL2

	readcol, '/Volumes/Data2/orelse/lemaux/deimos/N5281/cluster_catalog/NIRSPEC/starcats/setup' + setup + '/NIRSPECtargetsnoEWcat.cat', format='D,A,A,D,D,D,D,D', IDtar, mask, slit, RAtar, dectar, rbandtar, ibandtar, zbandtar	;get all target properties

	N = n_elements(RAtar)

	if N lt 3 then RAtartemp=dblarr(3)	;making RAtar a 3x1 arrayy to determine how many targets there are 
	if N lt 3 then for i=0, 3-1 do begin
	if i le N-1 then RAtartemp[i] = RAtar[i] else RAtartemp[i] = RAtartemp[i]
	endfor
	if N lt 3 then RAtar = RAtartemp

	readcol, '/Volumes/Data2/orelse/lemaux/deimos/N5281/cluster_catalog/NIRSPEC/potential.targets.NIRSPEC.cat', format='A,A,D,D', compmask1, compmask2, tartard, tartarPA	;get all relative positions and PAs from potential target list

	Ncomp = n_elements(compmask)
	
	if N gt 1 then index = intarr(N-1) else index = intarr(1)
	if N gt 1 then tartarPAforslit = dblarr(N-1) else tartarPAforslit = dblarr(1) 
	if N gt 1 then tartardforslit = dblarr(N-1) else tartardforslit = dblarr(1)

	checktar = where(RAtar gt 0)
	if n_elements(checktar) gt 1 then begin
	for i=1, N-1 do begin 	;get the PAs from the potential target list for actual targets

		index[i-1] = where(mask[0]+ '.' + slit[0] eq compmask1 and mask[i]+ '.' + slit[i] eq compmask2)
		tartardforslit[i-1] = tartard[index[i-1]]
		tartarPAforslit[i-1] = tartarPA[index[i-1]]*2.*!PI/360 
	endfor	
	endif
	
	if n_elements(checktar) eq 1 then begin ;if there is  only one target set the PA to be the galaxy's PA from the LFC cat 
	readcol, '/Volumes/Data2/orelse/lemaux/deimos/N5281/cluster_catalog/NIRSPEC/spectroscopic.roy.nep5281.clustermembers.cat', format='A,A,A', specID, specmask, specslit
	readcol, '/Volumes/Data2/orelse/lemaux/deimos/sc1604/cluster_catalog/EWcats/newEWcats/NIRSPEC/final.lfcpluscosmic.catalog.IDswPAs.cat', format='A,D', refID, SEXPA
	indexspec = where(mask[0] eq specmask and slit[0] eq specslit)
	tarspecID = specID[indexspec]
	indexforsex = where(refID eq tarspecID[0])
	tartarPAforslit[0] = SEXPA[indexforsex]
	tartardforslit[0] = 0.
	if tartarPAforslit[0] lt 0 then tartarPAforslit[0] = 360+tartarPAforslit[0] else tartarPAforslit[0] = tartarPAforslit[0]
	tartarPAforslit[0] = tartarPAforslit[0]*2*!PI/360
	endif	


	refslit = strarr(100)
	LFCID = strarr(100)
	dgs = dblarr(100)
	PAgs = dblarr(100)
	rbandgs = dblarr(100)
	ibandgs = dblarr(100)
	zbandgs = dblarr(100)

	for i=0, N-1 do begin	;get the PA and relative positions for the all guide stars for the setup calculated by starfind.pro
	readcol, '/Volumes/Data2/orelse/lemaux/deimos/N5281/cluster_catalog/NIRSPEC/starcats/setup' + setup + '/' + mask[i] + '.' + slit[i] + '.staroutput.cat', format = 'A,A,D,D,D,D,D', refslitindv, LFCIDindv, dgsindv, PAgsindv, rbandgsindv, ibandgsindv, zbandgsindv
	refslit = [refslit,refslitindv]
	LFCID = [LFCID, LFCIDindv]
	dgs = [dgs,dgsindv]
	PAgs = [PAgs,PAgsindv]
	rbandgs = [rbandgs,rbandgsindv]
	ibandgs = [ibandgs,ibandgsindv]
	zbandgs = [zbandgs,zbandgsindv]

	endfor

	Ngs = n_elements(PAgs)	

	


	if N ne 1 then RAcenter = (RAtar[0] + RAtar[1])/2 else RAcenter = RAtar[0]	;get the RA/dec center of the slit
	if N ne 1 then deccenter = (dectar[0] + dectar[1])/2 else deccenter = dectar[0]
	
	cdmatrix = dblarr(2,2)
	cdmatrix[0,0] = aa
	cdmatrix[0,1] = ab
	cdmatrix[1,0] = ba
	cdmatrix[1,1] = bb
	cdinvmatrix = invert(cdmatrix)
	
	;full spherical calculation of PA, d, d_RA, d_dec from intial point from fits header to center of slit
	p = atan( (-cos(deccenter*2*!PI/360)*sin( (CRVAL1-RAcenter)*2*!PI/360 ))/(cos(CRVAL2*2*!PI/360)*sin(deccenter*2*!PI/360)-sin(CRVAL2*2*!PI/360)*cos(deccenter*2*!PI/360)*cos( 2*!PI/360*(CRVAL1 - RAcenter))))
	if p lt 0 and CRVAL1 gt RAcenter then p = 2*!PI + p else if p lt 0 and CRVAL1 lt RAcenter then p = !PI + p else if p gt 0 and CRVAL1 gt RAcenter then p = !PI + p else if p gt 0 and CRVAL1 lt RAcenter then p = p

	deltadis = 360/(2*!PI)*acos(sin(deccenter*2*!PI/360)*sin(CRVAL2*2*!PI/360)+cos(CRVAL2*2*!PI/360)*cos(deccenter*2*!PI/360)*cos(2*!PI/360*(CRVAL1 - RAcenter)))
	deltaRA = deltadis*sin(p)
	deltadec = deltadis*cos(p)

	coormatrix =dblarr(2)
	coormatrix[0] = deltaRA
	coormatrix[1] = deltadec	

	;transform between difference in Ra/Dec to x/y pixels
	newpixarray = cdinvmatrix # coormatrix

	newxcenter = floor(refxpix + newpixarray[0])
	newycenter = floor(refypix + newpixarray[1])
	
	;1000x1000 guide image can change value if needed, corresponds to 180"x180" 
	EWedge = [newxcenter-500., newxcenter+500.]
	NSedge = [newycenter-500., newycenter+500.]

	readcol, '/Volumes/Data2/orelse/lemaux/deimos/N5281/cluster_catalog/NIRSPEC/starcats/setup'+ setup + '/goodguidestars.cat', format='A', gsID 	;get LFC IDs of good guide stars from list
	readcol, '/Volumes/Data2/orelse/lemaux/deimos/N5281/cluster_catalog/NIRSPEC/nep5281.lfc.IDradecmag.cat', format= 'A,D,D', LFCcatID, LFCRA, LFCdec 	;get RA/Dec of all LFC+COS objects from final LFC cat
	
	NLFC = n_elements(LFCcatID)
	
	gsRA = dblarr(3)
	gsdec = dblarr(3)
	PAgoodgs = dblarr(3)
	dgoodgs = dblarr(3)
	rbandgoodgs = dblarr(3)
	ibandgoodgs = dblarr(3)
	zbandgoodgs = dblarr(3)
	goodgsID = strarr(3)
	rank = indgen(3)+1 
	
	if RAtar[1] ne 0 then begin
	PAgoodgs2 = dblarr(3)
	dgoodgs2 = dblarr(3)
	endif else PAgoodgs2=dblarr(1)

	if RAtar[2] ne 0 then begin
	PAgoodgs3 = dblarr(3)
	dgoodgs3 = dblarr(3)	
	endif else PAgoodgs3=dblarr(1)

	for i=0, 2 do begin	;get the RA/Dec of the good guide stars and then get PA/d relative to the 1st target (should change to both or arbitrary #) and mag info

		for j=0, NLFC-1 do begin

		if LFCcatID[j] eq gsID[i] then gsRA[i] = LFCRA[j]
		if LFCcatID[j] eq gsID[i] then gsdec[i] = LFCdec[j]
		endfor 

		for j=0, Ngs-1 do begin
		if LFCID[j] eq gsID[i] and refslit[j] eq mask[0] + '.' + slit[0] then goodgsID[i] = LFCID[j]
		if LFCID[j] eq gsID[i] and refslit[j] eq mask[0] + '.' + slit[0] then PAgoodgs[i] = PAgs[j]
		if LFCID[j] eq gsID[i] and refslit[j] eq mask[0] + '.' + slit[0] then dgoodgs[i] = dgs[j]
		if LFCID[j] eq gsID[i] and refslit[j] eq mask[0] + '.' + slit[0] then rbandgoodgs[i] = rbandgs[j]
		if LFCID[j] eq gsID[i] and refslit[j] eq mask[0] + '.' + slit[0] then ibandgoodgs[i] = ibandgs[j]
		if LFCID[j] eq gsID[i] and refslit[j] eq mask[0] + '.' + slit[0] then zbandgoodgs[i] = zbandgs[j]
		endfor

		if n_elements(PAgoodgs2) gt 1 then for j=0, Ngs-1 do begin
		if LFCID[j] eq gsID[i] and refslit[j] eq mask[1] + '.' + slit[1] then PAgoodgs2[i] = PAgs[j]
		if LFCID[j] eq gsID[i] and refslit[j] eq mask[1] + '.' + slit[1] then dgoodgs2[i] = dgs[j]
		endfor

		if n_elements(PAgoodgs3) gt 1 then for j=0, Ngs-1 do begin
                if LFCID[j] eq gsID[i] and refslit[j] eq mask[2] + '.' + slit[2] then PAgoodgs3[i] = PAgs[j]
                if LFCID[j] eq gsID[i] and refslit[j] eq mask[2] + '.' + slit[2] then dgoodgs3[i] = dgs[j]
                endfor	
	endfor  
	
	dRAgoodgs = dgoodgs*sin(PAgoodgs*2*!PI/360)
	ddecgoodgs = dgoodgs*cos(PAgoodgs*2*!PI/360)
	if n_elements(PAgoodgs2) gt 1 then dRAgoodgs2 = dgoodgs2*sin(PAgoodgs2*2*!PI/360)	
	if n_elements(PAgoodgs2) gt 1 then ddecgoodgs2 = dgoodgs2*cos(PAgoodgs2*2*!PI/360)
	if n_elements(PAgoodgs3) gt 1 then dRAgoodgs3 = dgoodgs3*sin(PAgoodgs3*2*!PI/360)
        if n_elements(PAgoodgs3) gt 1 then ddecgoodgs3 = dgoodgs3*cos(PAgoodgs3*2*!PI/360)

	goodgsgsd = dblarr(3)
	goodgsgsPA = dblarr(3)
	goodgsgsdRA = dblarr(3)
	goodgsgsddec = dblarr(3)
	for i=0, 1 do begin ;calculate offsets of guide stars relative to one another
	
	goodgsgsd[i] = 360/(2*!PI)*acos(sin(gsdec[0]*2*!PI/360)*sin(gsdec[i+1]*2*!PI/360)+cos(gsdec[0]*2*!PI/360)*cos(gsdec[i+1]*2*!PI/360)*cos(2*!PI/360*( gsRA[0]- gsRA[i+1])))
	goodgsgsPA[i] = atan( (-cos(gsdec[i+1]*2*!PI/360)*sin( (gsRA[0]-gsRA[i+1])*2*!PI/360 ))/(cos(gsdec[0]*2*!PI/360)*sin(gsdec[i+1]*2*!PI/360)-sin(gsdec[0]*2*!PI/360)*cos(gsdec[i+1]*2*!PI/360)*cos( 2*!PI/360*(gsRA[0] -gsRA[i+1]))))	
	if goodgsgsPA[i] lt 0 and gsRA[0] gt gsRA[i+1] then goodgsgsPA[i] = 2*!PI + goodgsgsPA[i] else if goodgsgsPA[i] lt 0 and gsRA[0] lt gsRA[i+1] then goodgsgsPA[i] = !PI + goodgsgsPA[i] else if goodgsgsPA[i] gt 0 and gsRA[0] gt gsRA[i+1] then goodgsgsPA[i] = !PI + goodgsgsPA[i] else if goodgsgsPA[i] gt 0 and gsRA[0] lt gsRA[i+1] then goodgsgsPA[i] = goodgsgsPA[i]
	
	endfor 
	
	goodgsgsd[2] = 360/(2*!PI)*acos(sin(gsdec[1]*2*!PI/360)*sin(gsdec[2]*2*!PI/360)+cos(gsdec[1]*2*!PI/360)*cos(gsdec[2]*2*!PI/360)*cos(2*!PI/360*( gsRA[1]- gsRA[2])))
	goodgsgsPA[2] = atan( (-cos(gsdec[2]*2*!PI/360)*sin( (gsRA[1]-gsRA[2])*2*!PI/360 ))/(cos(gsdec[1]*2*!PI/360)*sin(gsdec[2]*2*!PI/360)-sin(gsdec[1]*2*!PI/360)*cos(gsdec[2]*2*!PI/360)*cos( 2*!PI/360*(gsRA[1] -gsRA[2]))))
	if goodgsgsPA[2] lt 0 and gsRA[1] gt gsRA[2] then goodgsgsPA[2] = 2*!PI + goodgsgsPA[2] else if goodgsgsPA[2] lt 0 and gsRA[1] lt gsRA[2] then goodgsgsPA[2] = !PI + goodgsgsPA[2] else if goodgsgsPA[2] gt 0 and gsRA[1] gt gsRA[2] then goodgsgsPA[2] = !PI + goodgsgsPA[2] else if goodgsgsPA[2] gt 0 and gsRA[1] lt gsRA[2] then goodgsgsPA[2] = goodgsgsPA[2]
		
	for i=0, 3-1 do begin 

	goodgsgsdRA[i] = 3600*goodgsgsd[i]*sin(goodgsgsPA[i])
	goodgsgsddec[i] = 3600*goodgsgsd[i]*cos(goodgsgsPA[i])
	endfor

	pgs = dblarr(3)
	deltadisgs = dblarr(3)
	deltarags = dblarr(3)
	deltadecgs = dblarr(3)
	dummy = dblarr(2,3)	

	for i=0, 2 do begin	;calculate PA and d for good guide stars from center of image
	pgs[i] = atan( (-cos(gsdec[i]*2*!PI/360)*sin( (RAcenter-gsRA[i])*2*!PI/360 ))/(cos(deccenter*2*!PI/360)*sin(gsdec[i]*2*!PI/360)-sin(deccenter*2*!PI/360)*cos(gsdec[i]*2*!PI/360)*cos( 2*!PI/360*(RAcenter - gsRA[i]))))
	deltadisgs[i] = 360/(2*!PI)*acos(sin(gsdec[i]*2*!PI/360)*sin(deccenter*2*!PI/360)+cos(deccenter*2*!PI/360)*cos(gsdec[i]*2*!PI/360)*cos(2*!PI/360*(RAcenter - gsRA[i])))

	if pgs[i] lt 0 and RAcenter gt gsRA[i] then pgs[i] = 2*!PI + pgs[i] else if pgs[i] lt 0 and RAcenter lt gsRA[i] then pgs[i] = !PI + pgs[i] else if pgs[i] gt 0 and RAcenter gt gsRA[i] then pgs[i] = !PI + pgs[i] else if pgs[i] gt 0 and RAcenter lt gsRA[i] then pgs[i] = pgs[i]

	deltarags[i] = deltadisgs[i]*sin(pgs[i])
	deltadecgs[i] =  deltadisgs[i]*cos(pgs[i])

	dummy[0:1,i] = [deltarags[i], deltadecgs[i]]
	endfor

	hmm = cdinvmatrix ## transpose(dummy)	;change from RA/dec ifference to pixel space 
	gspix = hmm + 500.

	
	ptar = dblarr(N)
	deltadistar = dblarr(N)
	deltaratar = dblarr(N)
	deltadectar = dblarr(N)
	dummytar = dblarr(2,N)

	if N eq 1 then begin
	dummytar[0:1,0] = [0, 0]
	hmmtar = cdinvmatrix ## transpose(dummytar)
	tarpix = hmmtar+500.
	endif 

	if N ne 1 then begin 
	
	for i=0, N-1 do begin ;get PA/d for target galaxies relative to the center of the image
	ptar[i] = atan( (-cos(dectar[i]*2*!PI/360)*sin( (RAcenter-RAtar[i])*2*!PI/360 ))/(cos(deccenter*2*!PI/360)*sin(dectar[i]*2*!PI/360)-sin(deccenter*2*!PI/360)*cos(dectar[i]*2*!PI/360)*cos( 2*!PI/360*(RAcenter - RAtar[i]))))
        deltadistar[i] = 360/(2*!PI)*acos(sin(dectar[i]*2*!PI/360)*sin(deccenter*2*!PI/360)+cos(deccenter*2*!PI/360)*cos(dectar[i]*2*!PI/360)*cos(2*!PI/360*(RAcenter - RAtar[i])))
	
	if ptar[i] lt 0 and RAcenter gt RAtar[i] then ptar[i] = 2*!PI + ptar[i] else if ptar[i] lt 0 and RAcenter lt RAtar[i] then ptar[i] = !PI + ptar[i] else if ptar[i] gt 0 and RAcenter gt RAtar[i] then ptar[i] = !PI + ptar[i] else if ptar[i] gt 0 and RAcenter lt RAtar[i] then ptar[i] = ptar[i]

	

        deltaratar[i] = deltadistar[i]*sin(ptar[i])
        deltadectar[i] =  deltadistar[i]*cos(ptar[i])

        dummytar[0:1,i] = [deltaratar[i], deltadectar[i]]
	endfor

	hmmtar = cdinvmatrix ## transpose(dummytar)	;transform from Ra/dec difference to pixel space
	tarpix = hmmtar+500.	
	endif

	file = '/Volumes/Data2/orelse/lemaux/deimos/N5281/cluster_catalog/NIRSPEC/starcats/setup' + setup + '/cutout.setup' + setup + '.fits' ;read in the finder chart
	finder = mrdfits(file,0)

	if n_elements(finder) lt 10 then begin ;if finder chart doesn't exist then make one

        mm = mrdfits('/Volumes/Data2/orelse/lemaux/LFC/N5281/nep5281_deep_norm.fits')

        finder = mm(EWedge[0]:EWedge[1], NSedge[0]:NSedge[1])
	
	mwrfits, finder, file

	finder = mrdfits(file,0)
	endif

	;readcol, '/Volumes/Data2/orelse/lemaux/deimos/sc1604/cluster_catalog/EWcats/newEWcats/NIRSPEC/starcats/' + setup + '/' + setup + 'description.cat', format='A', check 	

	;if n_elements(check) lt 2 then begin
	;print, n_elements(check)

	RAtarsexh = floor(RAtar/360*24) 
	RAtarsexm = floor((RAtar/360*24-floor(RAtar/360*24))*60)	
	RAtarsexs =  60*(((RAtar/360*24-floor(RAtar/360*24))*60)-floor(60*(RAtar/360*24-floor(RAtar/360*24) ) ))
	
	dectarsexdeg = floor(dectar)
	dectarsexam = floor((dectar-floor(dectar))*60)
	dectarsexas = 60*(((dectar-floor(dectar))*60)- floor(60*(dectar-floor(dectar))) )
	
	gsRAsexh = floor(gsRA/360*24)
        gsRAsexm = floor((gsRA/360*24-floor(gsRA/360*24))*60)
        gsRAsexs =  60*(((gsRA/360*24-floor(gsRA/360*24))*60)-floor(60*(gsRA/360*24-floor(gsRA/360*24) ) ))

        gsdecsexdeg = floor(gsdec)
        gsdecsexam = floor((gsdec-floor(gsdec))*60)
        gsdecsexas = 60*(((gsdec-floor(gsdec))*60)- floor(60*(gsdec-floor(gsdec))) )	

	tab = '09'XB
	
	; all this section is just to print out an ascii of the setup target and GS properties
	if keyword_set(maketable) then begin
	openw, lun, '/Volumes/Data2/orelse/lemaux/deimos/N5281/cluster_catalog/NIRSPEC/starcats/setup' + setup + '/description.setup' + setup + '.cat', /get_lun

	printf, lun, '#### Setup' + setup + ' description####', format='(a30)'
	printf, lun, '   ', format='(a8)'
	printf, lun, '####Target Info####', format = '(a20)'
	printf, lun, 'PA Tar-Tar',tab, 'd Tar-Tar(")',tab, 'd_RA Tar-Tar(")',tab, 'd_dec Tar-Tar(")', format='(a10, a1, a13, a1, a15, a1, a17)'
	printf, lun, tartarPAforslit[0]*360/(2*!PI), tab, tartardforslit[0], tab, tartardforslit[0]*sin(tartarPAforslit[0]), tab, tartardforslit[0]*cos(tartarPAforslit[0]), format = '(d7.3, a1, d7.3, a1, d7.3, a1, d7.3)'
	printf, lun, '   ', format='(a8)'	
	printf, lun, 'Name', tab, 'RA', tab, 'Dec', tab, 'r band', tab, 'i band', tab, 'z band', format='(a8, a1, a6, a1, a7, a1, a6, a1, a6, a1, a6)' 
		
	for i=0, N-1 do begin
	printf, lun, mask[i]+"."+slit[i], tab, RAtarsexh[i], ':', RAtarsexm[i], ':', RAtarsexs[i], tab, dectarsexdeg[i], ':', dectarsexam[i], ':', dectarsexas[i], tab, rbandtar[i], tab, ibandtar[i], tab, zbandtar[i], format='(a12, a1, i2.2, a1, i2.2, a1, d6.3, a1, i2.2, a1, i2.2, a1, d5.2, a1, d7.3, a1, d7.3, a1, d7.3)'
	endfor

	printf, lun, '    ', format='(a8)'
        printf, lun, '####GS-Target Movements####', format='(a30)'
	printf, lun, '    ', format='(a8)'
	printf, lun, 'GS Name', tab, 'LFC Name', tab, 'Rank', tab, 'RefSlit', tab, 'RA', tab, 'Dec', tab, 'r band', tab, 'd GS-Tar(")', tab, 'PA GS-Tar', tab, 'd_RA(")', tab, 'd_dec(")', format='(a10, a1, a10, a1, a10, a1, a10, a1, a10, a1, a10, a1, a10, a1, a12, a1, a12, a1, a12, a1, a12)'

	gsnum = strarr(3)
	gsnum[0] = 1 
	gsnum[1] = 2
	gsnum[2] = 3
	
	for i=0, 3-1 do begin			
	printf, lun, '52GS'+gsnum[i]+'NIR'+ setup, tab, goodgsID[i], tab, rank[i], tab, mask[0]+"."+slit[0], tab,  gsRAsexh[i], ':', gsRAsexm[i], ':', gsRAsexs[i], tab, gsdecsexdeg[i], ':', gsdecsexam[i], ':', gsdecsexas[i], tab, rbandgoodgs[i], tab, dgoodgs[i], tab, PAgoodgs[i], tab, dRAgoodgs[i], tab, ddecgoodgs[i], format='(a17, a1, a15, a1, i1, a1, a12, a1, i2.2, a1, i2.2, a1, d6.3, a1, i2.2, a1, i2.2, a1, d5.2, a1, d7.3, a1, d6.3, a1, d7.3, a1, d7.3, a1, d7.3)' 
	endfor

	if n_elements(PAgoodgs2) gt 1 then for i=0, 3-1 do begin
	printf, lun, '52GS'+gsnum[i]+'NIR'+ setup, tab, goodgsID[i], tab, rank[i], tab, mask[1]+"."+slit[1], tab, gsRAsexh[i], ':', gsRAsexm[i], ':', gsRAsexs[i], tab, gsdecsexdeg[i], ':', gsdecsexam[i], ':', gsdecsexas[i], tab, rbandgoodgs[i], tab, dgoodgs2[i], tab, PAgoodgs2[i], tab, dRAgoodgs2[i], tab, ddecgoodgs2[i], format='(a17, a1, a15, a1, i1, a1, a12, a1, i2.2, a1, i2.2, a1, d6.3, a1, i2.2, a1, i2.2, a1, d5.2, a1, d7.3, a1, d6.3, a1, d7.3, a1, d7.3, a1, d7.3)' 
	endfor	
	if n_elements(PAgoodgs3) gt 1 then for i=0, 3-1 do begin
        printf, lun, '52GS'+gsnum[i]+'NIR'+ setup, tab, goodgsID[i], tab, rank[i], tab, mask[2]+"."+slit[2], tab, gsRAsexh[i], ':', gsRAsexm[i], ':', gsRAsexs[i], tab, gsdecsexdeg[i], ':', gsdecsexam[i], ':', gsdecsexas[i], tab, rbandgoodgs[i], tab, dgoodgs3[i], tab, PAgoodgs3[i], tab, dRAgoodgs3[i], tab, ddecgoodgs3[i], format='(a17, a1, a15, a1, i1, a1, a12, a1, i2.2, a1, i2.2, a1, d6.3, a1, i2.2, a1, i2.2, a1, d5.2, a1, d7.3, a1, d6.3, a1, d7.3, a1, d7.3, a1, d7.3)' 
        endfor
	
	printf, lun, '    ', format='(a8)'
	printf, lun, '####GS-GS Movements####', format='(a25)'	
	printf, lun, 'Movement', tab, 'Ref GS Name', tab, 'To GS Name', tab, 'd GS-GS(")', tab, 'PA GS-GS', tab, 'd_RA(")', tab, 'd_dec(")', format='(a10, a1, a11, a1, a10, a1, a10, a1, a10, a1, a10, a1, a10, a1, a12, a1, a12, a1, a12, a1, a12)'

	printf, lun, gsnum[0] + '-' + gsnum[1], tab, '52GS'+gsnum[0]+'NIR'+ setup, tab, '52GS'+gsnum[1]+'NIR' + setup, tab, 3600*goodgsgsd[0], tab, 360/(2*!PI)*goodgsgsPA[0], tab, goodgsgsdRA[0], tab, goodgsgsddec[0], format='(a20, a1, a17, a1, a17, a1, d7.3, a1, d7.3, a1, d7.3, a1, d7.2)'	
	printf, lun, gsnum[0] + '-' + gsnum[2], tab, '52GS'+gsnum[0]+'NIR'+ setup, tab, '52GS'+gsnum[2]+'NIR'+ setup, tab, 3600*goodgsgsd[1], tab, 360/(2*!PI)*goodgsgsPA[1], tab, goodgsgsdRA[1], tab, goodgsgsddec[1], format='(a20, a1, a17, a1, a17, a1, d7.3, a1, d7.3, a1, d7.3, a1, d7.2)'
	printf, lun, gsnum[1] + '-' + gsnum[2], tab, '52GS'+gsnum[1]+'NIR'+ setup, tab, '52GS'+gsnum[2]+'NIR'+ setup, tab, 3600*goodgsgsd[2], tab, 360/(2*!PI)*goodgsgsPA[2], tab, goodgsgsdRA[2], tab, goodgsgsddec[2], format='(a20, a1, a17, a1, a17, a1, d7.3, a1, d7.3, a1, d7.3, a1, d7.2)'
	
	free_lun, lun
	endif

	finder = finder*(-1)
	if keyword_set(viewchart) then begin ; to display and make a hardcopy of the finder chart
	window, /free, xsize=1000, ysize=1000
	device, decomposed=1
	help, finder, /struct
	tv, finder	;display the finder chart

	for i=0, 2 do begin	;for all good guide stars draw a circle of radius 2.7" at their location

	tvcircle, 15., gspix[i,0], gspix[i,1], color='FFFF00'XL, THICK=3 
	
	endfor

	for i=0, N-1 do begin	;for all targets draw a circle of radius 1.8" at their location
	
	tvcircle, 10., tarpix[i,0.], tarpix[i,1.], color='FF0000'XL, THICK=2

	endfor		

	slitsize = dblarr(2)	;2 array with slit width and slit length 1" x 42"
	slitsize[0] = 1./0.18
	slitsize[1] = 42./0.18	
	
	tvbox, slitsize, 499., 499., color='0000FF'XL, angle=-tartarPAforslit[0]/(2*!PI)*360	;draw a slit at the center of the image with PA given by the PA calculated between targets
	device, decomposed=0
	endif

	if keyword_set(makechart) then begin	
	entry_device = !d.name
	set_plot, 'PS'
	device, /color, bits_per_pixel=24, filename = '/Volumes/Data2/orelse/lemaux/deimos/N5281/cluster_catalog/NIRSPEC/starcats/setup' + setup + '/setup'+ setup+ '.finder.ps', xoffset=0.75, yoffset= 0.75, xsize=6.5, ysize=9, /inches
	
	loadcolors			;load custom color table with loadcolors.pro
	windowXSize = !d.x_size		; get the x-window size in device coordinates of the PS window
   	windowYSize = !d.y_size		
   	x0 = (0.5/6.5) * windowXSize	; starting x values is 0.5 in inches traqnsform to device coordinates
   	y0 = (4.25/9.0) * windowYSize
   	xsize = (5.75/6.5) * windowXSize	;the image is 5.75 inches inside the 6.5 inch window, transform to device coordinates
   	ysize = (5.75/9.0) * windowYSize
	tv, bytscl(finder, min=-13.3, max=-11.7), x0, y0, XSIZE=xsize, YSIZE=ysize 

	for i=0, 2 do begin     ;for all good guide stars draw a circle of radius 2.7" at their location

        tvcircle, 15.*xsize/1000, gspix[i,0]*xsize/1000+x0, gspix[i,1]*ysize/1000+y0, color=6, THICK=3	;if I change the pixel size of the original cutout then 1000 needs to be changed in here

        endfor

        for i=0, N-1 do begin   ;for all targets draw a circle of radius 1.8" at their location

        tvcircle, 10.*xsize/1000, tarpix[i,0.]*xsize/1000+x0, tarpix[i,1.]*ysize/1000+y0, color=5, THICK=2

        endfor

        slitsize = dblarr(2)    ;2 array with slit width and slit length 1" x 42"
        slitsize[0] = 1./0.18*xsize/1000
        slitsize[1] = 42./0.18*ysize/1000	;pixel scale is 0.18 in LFC unbinned differed if using a different image

        tvbox, slitsize, 499.*xsize/1000+x0, 499.*ysize/1000+y0, color=4, angle=-tartarPAforslit[0]/(2*!PI)*360, THICK=2	;this will also need to be changed if original image size is change 499 to new center

	for i=0, 2 do begin
	xyouts, gspix[i,0]*xsize/1000-x0, gspix[i,1]*ysize/1000+y0+300, rank[i], charsize=1.5, color=6, /device		;label the guide stars with their priority on the image
	endfor	

	for i=0, N-1 do begin
	xyouts, tarpix[i,0.]*xsize/1000+x0+300, tarpix[i,1.]*ysize/1000+y0, mask[i] + '.' + slit[i], charsize=1, color=5, /device 	;label the target names
	endfor

	linesize = dblarr(2)			;for creating N E orientation on image
	linesize[0] = .1/0.18*xsize/1000
	linesize[1] = 15./0.18*ysize/1000	
	tvbox, linesize, 200.*xsize/1000+x0, 900.*ysize/1000+y0, angle=0, THICK=2, color=6
	tvbox, linesize, 200.*xsize/1000+x0-linesize[1]/2, 900.*ysize/1000+y0-linesize[1]/2, color=6, angle=90, THICK=3
	xyouts, 200.*xsize/1000+x0 - 150,  900.*ysize/1000+y0 + linesize[1]/2 + 200, 'N', charsize=1.25, color=6, /device
	xyouts, 200.*xsize/1000+x0 - linesize[1] - 300,  900.*ysize/1000+y0 - linesize[1]/2-150, 'E', charsize=1.25, color=6, /device

	xyouts, (0.5/7.0)*windowXsize, 10000, 'Setup ' + setup, charsize=1.75, color=0, /device 
	xyouts, (0.5/7.0)*windowXsize, 9500, 'Base Target is ' + mask[0] + '.' + slit[0], charsize=1, color=0, /device
	xyouts, (0.15/7.0)*windowXsize, 8500, 'PA Target-Target   d Target-Target(")   d_RA(")   d_dec(")', charsize=1, color=0, /device
	xyouts, (0.15/7.0)*windowXsize, 8000, string(tartarPAforslit[0]*360/(2*!PI), format='(D6.2)') + '              ' + string(tartardforslit[0], format='(D6.3)') + '              ' + string(tartardforslit[0]*sin(tartarPAforslit[0]), format='(D7.3)') + '   ' + string(tartardforslit[0]*cos(tartarPAforslit[0]), format='(D7.3)'), charsize=1, color=0, /device	

	xyouts, (0.15/7.0)*windowXsize, 7000, 'Name         RA               Dec             r-band   i-band   z-band', charsize=1, color=0, /device

	indexarr = indgen(100)

	;for i=0, 0 do begin
	;xyouts, (0.15/7.0)*windowXsize, 6500-500*indexarr[i], mask[i]+'.'+slit[i] + '  ' +  string(RAtarsexh[i], format='(D4.0)') + ':' + string(RAtarsexm[i], format='(D3.0)') + ':0' + string(RAtarsexs[i], format='(D5.3)')  + '  ' + string(dectarsexdeg[i], format='(D4.0)') + ':' + string(dectarsexam[i], format='(D3.0)') + ':' + string(dectarsexas[i], format='(D6.3)') + '  ' + string(rbandtar[i], format='(D7.3)') + '  ' +  string(ibandtar[i], format='(D7.3)') + '  '  + string(zbandtar[i], format='(D7.3)'), color=0, charsize=1, /device
	;endfor

	for i=0, N-1 do begin
        xyouts, (0.15/7.0)*windowXsize, 6500-500*indexarr[i], mask[i]+'.'+slit[i] + '  ' +  string(RAtarsexh[i], format='(D4.0)') + ':' + string(RAtarsexm[i], format='(D3.0)') + ':' + string(RAtarsexs[i], format='(D6.3)')  + '  ' + string(dectarsexdeg[i], format='(D4.0)') + ':' + string(dectarsexam[i], format='(D3.0)') + ':' + string(dectarsexas[i], format='(D6.3)') + '  ' + string(rbandtar[i], format='(D7.3)') + '  ' +  string(ibandtar[i], format='(D7.3)') + '  '  + string(zbandtar[i], format='(D7.3)'), color=0, charsize=1, /device
        endfor

	xyouts, (0.15/7.0)*windowXsize, 6500-500*indexarr[N-1+2], 'LFC Name  Rank  RefSlit       RA             Dec              r-band    d_RA(")   d_dec(")', color=0, charsize=1, /device

	for i=0, 3-1 do begin
	xyouts, (0.15/7.0)*windowXsize, 6500-500*indexarr[N-1+3+i], goodgsID[i] + '        ' + string(rank[i], format='(I1)') + '  ' + mask[0]+'.'+slit[0] + '  ' + string(gsRAsexh[i], format='(D4.0)') +':' + string(gsRAsexm[i], format='(D3.0)') + ':' + string(gsRAsexs[i], format='(D6.3)') + '  ' + string(gsdecsexdeg[i], format='(D4.0)') + ':' + string(gsdecsexam[i], format='(D3.0)') + ':' + string(gsdecsexas[i], format='(D6.3)') + '   ' + string(rbandgoodgs[i], format='(D7.3)') + '   ' + string(dRAgoodgs[i], format='(D7.3)') + '   ' + string(ddecgoodgs[i], format='(D7.3)'), color=0, charsize=1,  /device
	endfor

	;for i=2, 2 do begin
        ;xyouts, (0.15/7.0)*windowXsize, 6500-500*indexarr[N-1+3+i], goodgsID[i] + '        ' + string(rank[i], format='(I1)') + '  ' + mask[0]+'.'+slit[0] + '  ' + string(gsRAsexh[i], format='(D4.0)') +':' + string(gsRAsexm[i], format='(D3.0)') + ':' + string(gsRAsexs[i], format='(D6.3)') + '  ' + string(gsdecsexdeg[i], format='(D4.0)') + ':' + string(gsdecsexam[i], format='(D3.0)') + ':0' + string(gsdecsexas[i], format='(D5.3)') + '   ' + string(rbandgoodgs[i], format='(D7.3)') + '   ' + string(dRAgoodgs[i], format='(D7.3)') + '   ' + string(ddecgoodgs[i], format='(D7.3)'), color=0, charsize=1,  /device
        ;endfor

	;for i=2, 2 do begin
        ;xyouts, (0.15/7.0)*windowXsize, 6500-500*indexarr[N-1+3+i], goodgsID[i] + '  ' + string(rank[i], format='(I1)') + '  ' + mask[0]+'.'+slit[0] + '  ' + string(gsRAsexh[i], format='(D4.0)') +':0' + string(gsRAsexm[i], format='(D2.0)') + ':' + string(gsRAsexs[i], format='(D6.3)') + '  ' + string(gsdecsexdeg[i], format='(D4.0)') + ':' + string(gsdecsexam[i], format='(D3.0)') + ':0' + string(gsdecsexas[i], format='(D5.3)') + '   ' + string(rbandgoodgs[i], format='(D7.3)') + '   ' + string(dRAgoodgs[i], format='(D7.3)') + '   ' + string(ddecgoodgs[i], format='(D7.3)'), color=0, charsize=1,  /device
        ;endfor
	
	;for i=2, 2 do begin
	;xyouts, (0.15/7.0)*windowXsize, 6500-500*indexarr[N-1+3+i], goodgsID[i] + '  ' + string(rank[i], format='(I1)') + '  ' + mask[0]+'.'+slit[0] + '  ' + string(gsRAsexh[i], format='(D4.0)') +':0' + string(gsRAsexm[i], format='(D2.0)') + ':' + string(gsRAsexs[i], format='(D6.3)') + '  ' + string(gsdecsexdeg[i], format='(D4.0)') + ':' + string(gsdecsexam[i], format='(D3.0)') + ':0' + string(gsdecsexas[i], format='(D5.3)') + '   ' + string(rbandgoodgs[i], format='(D7.3)') + '   ' + string(dRAgoodgs[i], format='(D7.3)') + '   ' + string(ddecgoodgs[i], format='(D7.3)'), color=0, charsize=1,  /device
	;endfor	
	
	if n_elements(PAgoodgs2) gt 1 then N=N+3	;stopped after 2 targets, for three targets go to description file
	if n_elements(PAgoodgs2) gt 1 then for i=0, 3-1 do begin
	xyouts, (0.15/7.0)*windowXsize, 6500-500*indexarr[N-1+3+i], goodgsID[i] + '        ' + string(rank[i], format='(I1)') + '  ' + mask[1]+'.'+slit[1] + '  ' + string(gsRAsexh[i], format='(D4.0)') +':' + string(gsRAsexm[i], format='(D3.0)') + ':' + string(gsRAsexs[i], format='(D6.3)') + '  ' + string(gsdecsexdeg[i], format='(D4.0)') + ':' + string(gsdecsexam[i], format='(D3.0)') + ':' + string(gsdecsexas[i], format='(D6.3)') + '   ' + string(rbandgoodgs[i], format='(D7.3)') + '   ' + string(dRAgoodgs2[i], format='(D7.3)') + '  ' + string(ddecgoodgs2[i], format='(D7.3)'), color=0, charsize=1,  /device	
	endfor

	;if n_elements(PAgoodgs2) gt 1 then for i=2, 2 do begin
        ;xyouts, (0.15/7.0)*windowXsize, 6500-500*indexarr[N-1+3+i], goodgsID[i] + '        ' + string(rank[i], format='(I1)') + '  ' + mask[1]+'.'+slit[1] + '  ' + string(gsRAsexh[i], format='(D4.0)') +':' + string(gsRAsexm[i], format='(D3.0)') + ':' + string(gsRAsexs[i], format='(D6.3)') + '  ' + string(gsdecsexdeg[i], format='(D4.0)') + ':' + string(gsdecsexam[i], format='(D3.0)') + ':0' + string(gsdecsexas[i], format='(D5.3)') + '   ' + string(rbandgoodgs[i], format='(D7.3)') + '   ' + string(dRAgoodgs2[i], format='(D7.3)') + '  ' + string(ddecgoodgs2[i], format='(D7.3)'), color=0, charsize=1,  /device 
        ;endfor
	gsnum = indgen(3)+1	
	xyouts, (0.15/7.0)*windowXsize, 6500-500*[N-1+7], 'Movement   Ref GS Name   To GS Name   d_RA(")   d_dec(") ', charsize=1, color=0, /device	
	xyouts, (0.15/7.0)*windowXsize, 6500-500*[N-1+8], string(gsnum[0], format='(I1)') + '-' + string(gsnum[1], format='(I1)') + '       52GS'+ string(gsnum[0], format='(I1)') +'NIR'+ setup + '      52GS'+ string(gsnum[1], format='(I1)') +'NIR' + setup + '      ' +  string(goodgsgsdRA[0], format='(D8.3)') + '   ' +  string(goodgsgsddec[0], format='(D8.3)'), charsize=1, color=0, /device	
	xyouts, (0.15/7.0)*windowXsize, 6500-500*[N-1+9], string(gsnum[0], format='(I1)') + '-' + string(gsnum[2], format='(I1)') + '       52GS'+ string(gsnum[0], format='(I1)') +'NIR'+ setup + '      52GS'+ string(gsnum[2], format='(I1)') +'NIR' + setup + '      ' +  string(goodgsgsdRA[1], format='(D8.3)') + '   ' +  string(goodgsgsddec[1], format='(D8.3)'), charsize=1, color=0, /device
	xyouts, (0.15/7.0)*windowXsize, 6500-500*[N-1+10], string(gsnum[1], format='(I1)') + '-' + string(gsnum[2], format='(I1)') + '       52GS'+ string(gsnum[1], format='(I1)') +'NIR'+ setup + '      52GS'+ string(gsnum[2], format='(I1)') +'NIR' + setup + '      ' +  string(goodgsgsdRA[2], format='(D8.3)') + '   ' +  string(goodgsgsddec[2], format='(D8.3)'), charsize=1, color=0, /device

	device, /close_file	
	set_plot, entry_device
	endif
	
	

end
	 	
