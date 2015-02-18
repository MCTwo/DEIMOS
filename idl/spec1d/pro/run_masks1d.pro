;+
; NAME:
;
; RUN_MASKS1d
;
; PURPOSE:
;
;  launch the next n masks through spec1d
;
; CATEGORY:
;
;  spec1d control
;
; CALLING SEQUENCE:
;
;  run_masks1d [,n]
;
; INPUTS:
;
;   n -- number of masks (default: 5)
;
; RESTRICTIONS:
;
;  requires a file called spec1dcontrol.txt in ~/scripts containing
;   one or two lines: first, the version number to check for; and
;   second, the # of masks to run at once.  also required is a maskcontrol.txt
;   file containing 2 or 3 lines: first, the version number to check
;   for.
;    The location
;   of the spec2d scripts should be passed as the second line of this
;   file.  A number of masks can be listed as the third line of this
;   file, rather than passed at the command line.
;
; MODIFICATION HISTORY:
;   2003jul14 JAN
;-


pro run_masks1d,number,checkjobs=checkjobs

if n_elements(checkjobs) eq 0 then checkjobs=0


openr,2,'~/scripts/spec1dcontrol.txt'

version='abc'
n='1'
path='abc'
readf,2,version
readf,2,n
close,2


openr,2,'~/scripts/maskcontrol.txt'

version_2d='abc'
n_2d='1'
path='abc'
readf,2,version_2d
readf,2,n_2d
readf,2,path
close,2

if n_elements(number) eq 0 then $ 
   if n_elements(n) ne 0 then number=fix(n) else number=5

waitstring=''


; if desired, check to see if we might hold up any multiprocessor
; jobs; if so, submit in a fashion that waits for them to complete.
; The queue itself should make us wait for single-processor jobs.
; Also, in this case refrain from totally filling the queue.
if checkjobs then begin
    spawn,'qstat -R',qoutput
    qoutput=qoutput[5:*]
    njobs=n_elements(qoutput)
    jobnum=strarr(njobs)
    nnodes=intarr(njobs)
    status=strarr(njobs)
    for i=0,njobs-1 do begin
        split=strsplit(qoutput[i],/extract)
        jobnum[i]=(strsplit(split[0],'.',/extract))[0]
        nnodes[i]=split[3]
        status[i]=split[7]
    endfor
    holdup=nnodes gt 1 and (status eq 'Q' or status eq 'H')
    if total(holdup) gt 0 then begin
        wh=where(holdup,nhold)
        waitstring='-W depend=after'
        for i=0,nhold-1 do begin
            waitstring=waitstring+':'+jobnum[i]
        endfor                
        waitstring=waitstring+' '
    endif
    number=number < (20-njobs)

endif

path1d=strmid(path,0,strpos(path,'spec2d'))+'spec1d'


 search1d=concat_dir(path1d,'spec1d.*.200*.sh')


shfiles1d=findfile(search1d,count=shcount1d)

; make list of mask names in spec1d files
masks1d=strmid(findfile(search1d),strpos(search1d[0],'spec1d.') +7,99)

for i=0,n_elements(masks1d)-1 do begin  
    tmp=strsplit(masks1d[i],'.',/extract) 
    masks1d[i]=tmp[0] &$
endfor



donearr=intarr(shcount1d)
done1darr=intarr(shcount1d)
nowprocessing=intarr(shcount1d)
spawn,'qstat -f',qstat

for i=0,shcount1d-1 do begin

; check if this mask is currently running
    spawn,'grep "cd " '+shfiles1d[i],output
    spawn,'grep "PBS -N " '+shfiles1d[i],jobstring
    split=strsplit(output,/extract)
    if n_elements(split) ge 2 then maskpath=split[1] $
      else maskpath=''
    split=strsplit(jobstring,/extract)
      if n_elements(split) ge 3 then jobname=split[2] $
      else jobname=''
    
    if maskpath eq '' then begin
        donearr[i]=1 
        print,'No path found for script '+shfiles1d[i]
    endif else begin
        nowprocessing[i]=total(strpos(qstat,jobname) ge 0) ne 0 $
          AND jobname ne ''
        doneversion='abc'

; check if spec2d has been run here
        donefile=concat_dir(maskpath,'doneprocessing.txt')
        donefile2=findfile(donefile,count=doneexists)

        if doneexists ne 0 then begin
            openr,2,donefile2
            readf,2,doneversion
            close,2
            donearr[i]=doneversion eq version_2d
        endif

; check if spec1d has been run here
        donefile=concat_dir(maskpath,'donespec1d.txt')
        donefile2=findfile(donefile,count=doneexists)

        if doneexists ne 0 then begin
            openr,2,donefile2
            readf,2,doneversion
            close,2
            done1darr[i]=doneversion eq version
        endif

    endelse

endfor


     whnotdone=where(donearr eq 1 and done1darr eq 0 $
                     and nowprocessing eq 0,notdonect)


; submit masks
     if notdonect gt 0 then begin

         for i=0,(notdonect-1) < (number-1) do begin

         spawn,'qsub '+waitstring+shfiles1d[whnotdone[i]]
         print,'qsub '+waitstring+shfiles1d[whnotdone[i]]

         endfor
     endif
return
end


