; To make ascii catalog out of zspec .fits output files

pro tablefromzspec, user, mask   

	goo = 'zspec.' + user + '.' + mask + '.*.fits'
	foo = mrdfits(goo,1)
	openw, outfile, mask + '.objname.out', /get_lun, WIDTH = 120
	printf, outfile, foo.objname, FORMAT = '(1(A))'
	free_lun, outfile
	openw, outfile, mask + '.slitname.out', /get_lun, WIDTH = 120
        printf, outfile, foo.slitname, FORMAT = '(1(I))'
        free_lun, outfile
	openw, outfile, mask + '.maskname.out', /get_lun, WIDTH = 120
        printf, outfile, foo.maskname, FORMAT = '(1(A))'
        free_lun, outfile
	openw, outfile, mask + '.z.out', /get_lun, WIDTH = 120
        printf, outfile, foo.z, FORMAT = '(1(D))'
        free_lun, outfile
	openw, outfile, mask + '.z_err.out', /get_lun, WIDTH = 120
        printf, outfile, foo.z_err, FORMAT = '(1(D))'
        free_lun, outfile
	openw, outfile, mask + '.zquality.out', /get_lun, WIDTH = 120
        printf, outfile, foo.zquality, FORMAT = '(1(I))'
        free_lun, outfile
	openw, outfile, mask + '.comment.out', /get_lun, WIDTH = 120
        printf, outfile, foo.comment, FORMAT = '(1(A))'
        free_lun, outfile
end
