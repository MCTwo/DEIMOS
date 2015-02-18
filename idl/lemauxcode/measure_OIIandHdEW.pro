pro measure_OIIandHdEW, spec

	mm = mrdfits(spec,1)
	OII = measure_simple_ew(mm,$
	      [[3653, 3723], [3733.5,3803]],$
	      [3723,3733.5])
	Hd = measure_simple_ew(mm,$
	     [[4017.,4057.],[4153.,4193.]], $
	     [4083.5,4122.250])
	print, 'EW(OII) is:', OII, '    EW(Hd) is:', Hd

end 
