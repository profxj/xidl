;+
; NAME:
;
; RUN_MASKS
;
; PURPOSE:
;
;  launch the next n masks 
;
; CATEGORY:
;
;  spec2d control
;
; CALLING SEQUENCE:
;
;  run_masks [,n]
;
; INPUTS:
;
;   n -- number of masks (default: 5)
;
; RESTRICTIONS:
;
;  requires a file called maskcontrol.txt in ~/scripts containing
;   three lines: first, the version number to check for.  The location
;   of the spec2d scripts should be passed as the second line of this
;   file.  A number of masks can be listed as the third line of this
;   file, rather than passed at the command line.
;
; MODIFICATION HISTORY:
;   2003jul11 JAN
;-


pro run_masks,number


openr,2,'~/scripts/maskcontrol.txt'

version='abc'
n='1'
path='abc'
readf,2,version
readf,2,n
readf,2,path
close,2

if n_elements(number) eq 0 then $ 
   if n_elements(n) ne 0 then number=fix(n) else number=5

search=concat_dir(path,'spec2d.*.200*.sh')
shfiles=findfile(search,count=shcount)

donearr=intarr(shcount)
nowprocessing=intarr(shcount)
spawn,'qstat -f',qstat

for i=0,shcount-1 do begin

    spawn,'grep "cd " '+shfiles[i],output
    spawn,'grep "PBS -N " '+shfiles[i],jobstring
    split=strsplit(output,/extract)
    if n_elements(split) ge 2 then maskpath=split[1] $
      else maskpath=''
    split=strsplit(jobstring,/extract)
      if n_elements(split) ge 3 then jobname=split[2] $
      else jobname=''
    
    if maskpath eq '' then begin
        donearr[i]=1 
        print,'No path found for script '+shfiles[i]
    endif else begin
        nowprocessing[i]=total(strpos(qstat,jobname) ge 0) ne 0 $
          AND jobname ne ''
        doneversion='abc'
        donefile=concat_dir(maskpath,'doneprocessing.txt')
        donefile2=findfile(donefile,count=doneexists)
        spawn,'chmod g+w '+maskpath+' &',result,err
        spawn,'chmod g+w '+maskpath+'/* &',result,err
        spawn,'chmod g+w '+maskpath+'/*/* &',result,err

        if doneexists ne 0 then begin
            openr,2,donefile2
            readf,2,doneversion
            close,2
            donearr[i]=doneversion eq version
        endif

    endelse

endfor

     whnotdone=where(donearr eq 0 and nowprocessing eq 0,notdonect)
     if notdonect gt 0 then begin

         for i=0,(notdonect-1) < (number-1) do begin

         spawn,'qsub '+shfiles[whnotdone[i]]
         print,'qsub '+shfiles[whnotdone[i]]

         endfor
     endif
return
end


