;+
;; Script for running 2D co-adds on LRIS MOS data.
;;
;; This script should be '.run' in the directory level where the
;; wavelenght solutions and mask trace files reside, also a level below
;; 'Science' with the reduced data from the individual frames (from
;; running LONG_REDUCE as usual. 
;;
;; So far only tested for LRIS-B reductions
;;
;; KG Lee   06/07/2015 - Code adapted to run long_coadd2d_run and
;;                       long_reduce 
;-

cd, current=maskpath

maskpath = maskpath+'/'

check = 1
long_coadd2d_run, maskpath, check = check

planfil_orig = 'plan.par'
planfil_2d = 'plan_coadd2d.par'

if not file_test(planfil_orig) then begin
   read, planfil_orig, prompt='plan.par not found. Please enter name of plan file.'
endif

;; Read in original plan.par file. plan_coadd2d.par will be modified
;; from this
planstr = yanny_readone(planfil_orig, hdr=planhdr)

;; Modifying header in plan.par ...
;; Loop through header lines and change 'plan.log' to
;; 'plan_coadd2d.log' etc.
;;
;; Also change indir to local directory and strip comment lines with
;; '#' from planhdr

print, 'Modifying header of plan file'
nlines = n_elements(planhdr)
commentlines = []
for ii=0, nlines-1 do begin
   
   str1ch = strmid(strtrim(planhdr[ii],2),0,1)
   if strmatch(str1ch, '#') OR strmatch(strtrim(planhdr[ii],2),'') EQ 1 $
   then commentlines=[commentlines,ii]
   
   str1word = (strsplit(planhdr[ii], /extract))[0] 
   if strmatch(str1word, 'logfile') then $
      planhdr[ii] = "logfile  'plan_coadd2d.log' # Log file"
   if strmatch(str1word, 'plotfile') then $
      planhdr[ii] = "plotfile  'plan_coadd2d.ps' # Plot file"
   if strmatch(str1word, 'indir') then begin
      oldrawstr = (strsplit(planhdr[ii],/extract))[1]
      planhdr[ii] = "indir  './' # Raw data directory"
   endif
endfor
remove, commentlines, planhdr

;; Now modify input files, while checking that we have the files we
;; need in the local directory. In the case of the twiflat, we want to
;; copy it over from the Raw directory

;; Check for bias
;findbias = where(strmatch(planstr.flavor, 'bias'))
;biasfile = planstr[findbias].filename
;if not file_test('bias-'+biasfile) AND $
;   not file_test('superbias-'+biasfile) then begin
;   print, 'No superbias file found. Copying over single bias frame from raw'
;   findbias = where(strmatch(planstr.flavor, 'bias') EQ 1)
;   biasfile = planstr[findbias].filename
;   oldrawdir = (strsplit(oldrawstr, "'",/extract))[0]
;   print, biasfile
;   if not file_test(biasfile+'*') then file_copy, oldrawdir+biasfile+'*', '.'
;endif

;; Check for wavelength solution
print, 'Looking for wavelength solutions'
findarc = where(strmatch(planstr.flavor, 'arc'))
arcfile = planstr[findarc].filename
if strmatch(arcfile,'*.gz') then $
   arcfile=strmid(arcfile,0,strlen(arcfile)-3)
if not file_test('wave-'+arcfile) then begin
   if file_test('wave-'+arcfile+'.gz') then arcfile=arcfile+'.gz' else begin
      print, 'Error: No wavelength solution found'
      stop
   endelse
endif

;; Copy over twiflat
;findtwi = where(strmatch(planstr.flavor,'twiflat'))
;twifile = planstr[findtwi].filename
;oldrawdir = (strsplit(oldrawstr, "'",/extract))[0]
;stop
;if not file_test(twifile+'*') then file_copy, oldrawdir+twifile+'*' , '.'

;; Check for dome flat
print, 'Looking for dome flats'
finddome = where(strmatch(planstr.flavor,'domeflat'))
domefile = planstr[finddome].filename
if strmatch(domefile,'*.gz') then $
   domefile=strmid(domefile,0,strlen(domefile)-3)
if not file_test('pixflat-'+domefile) then begin
   if file_test('pixflat-'+domefile+'.gz') then domefile=domefile+'.gz' else begin
      print, 'Error: No pixflat file found'
      stop
   endelse
endif
   
print, 'Looking for science file'
;; Check for which sci file has the correspondingly named coadd
;; file. The rest get remooved
findsci = where(strmatch(planstr.flavor, 'sci*') EQ 1, nsci)

killsci = []
for ii=0, nsci-1 do begin
   filnametmp = planstr[findsci[ii]].filename
   if strmatch(filnametmp, '*.gz') then begin
      len=strlen(filnametmp)
      filnametmp = strmid(filnametmp, 0, len-3)
      print, filnametmp
   endif

   if file_test('coadd_2d-'+filnametmp) then $
      planstr[findsci[ii]].filename = $
      'coadd_2d-'+filnametmp else  killsci = [killsci, findsci[ii]]

endfor
remove, killsci, planstr

print, 'Writing new plan file'
yanny_write, planfil_2d, ptr_new(planstr), hdr=planhdr, /align

;; Run LONG_REDUCE with new plan file
long_reduce, planfil_2d

END
