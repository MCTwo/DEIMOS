pro get_lambda_central, spec, initlam, bandpass

	path = getenv('D2_RESULTS')
	mmB = mrdfits(path + 'sc1604/' + spec,3)
	mmR = mrdfits(path + 'sc1604/' + spec,4)
	if initlam gt mmB.lambda[4095] then mm = mmR else mm = mmB

	dental = where(mm.lambda gt initlam)
	plan = dental[0]
	lisa = max(mm.spec[(plan-bandpass/2):(plan+bandpass/2)])
	needs = where(mm.spec eq lisa)
	braces = mm.lambda[needs]

	if abs(braces-initlam) gt 5 then begin 
	message, "Initial wavelength guess too far from peak value, limit bandpass or change guess!!!", /info
	endif
	
 	if abs(braces-initlam) le 5 then begin
	print, braces
	endif
end 
