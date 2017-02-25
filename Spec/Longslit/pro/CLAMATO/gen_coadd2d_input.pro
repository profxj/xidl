;+
; Script to generate input file for COADD_CLAMATO_LRIS.PRO . Just go
; into slitmask input file, and create full list of objects ranked by
; descending order of y-position.
;
; (Designed) slits that fell in the CCD chip gap need to be indicated with the
; REM_SLIT keyword (as a vector of obj-id's), in order to match
; the slit number convention from long_reduce
;
; SLITS_FIL and OUTFIL are optional; default input and output files
; will be based on maskname.
;
; SCI_FIL is also optional, and is used to read in HDU5 of a
; sci*.fits.gz file in order to list the first extracted object IDs
; corresponding to each slit. 
;
; The output list will need to be edited manually in order to specify
; the actual traces we want to extract. 
;
; Example: gen_coadd2d_input, maskname, slits_fil=slits_fil, $
;                       outfil=outfil, rem_slit=rem_slit, sci_fil=sci_fil
;
;-
pro gen_coadd2d_input, maskname, slits_fil=slits_fil, $
                       outfil=outfil, rem_slit=rem_slit, $
                       sci_fil=sci_fil
  
  if not keyword_set(maskname) then begin
     print, 'Error: Please specify name of slitmask'
     stop
  endif
  
  if not keyword_set(outfil) then $
     outfil= 'coadd2d_input_'+maskname+'_Blue600.txt'

  if not keyword_set(slits_fil) then begin
     slits_fil= maskname+'.slits'
     slits_fil='~/lya/3d_recon/ts/clamato2016/slitmasks/'+slits_fil
  endif
     
  readcol, slits_fil, x_star, y_star, prior, name, mag, $
           f='f,f,x,x,x,x,x,l,a,f', skip=63

;; Remove alignment stars from list
starmatch = where(strmatch(name, 'star*') EQ 1)
remove, starmatch, x_star, y_star, prior, name, mag

;; Sort by reverse y-position
ysort = reverse(sort(y_star))
name  = name[ysort]
mag   = mag[ysort]
prior = prior[ysort]

;; *** Remove objects that did not get slits. Primarily slits that
;; fell in chip gap ***
if keyword_set(rem_slit) then begin
   trace_fail = rem_slit - 1
   remove, trace_fail, name, mag, prior 
endif

nslits = n_elements(name)
slitnum = findgen(nslits)+1

;; Obj ID defaults to slit number, but....
;; If a sci*.fits.gz file is defined, then read this to match the obj
;; IDs as extracted for each slit number. Just pick the first obj ID
;; per slit
objnum = slitnum
if keyword_set(sci_fil) then begin
   sciobj = mrdfits(sci_fil, 5)
   for ii=0, nslits-1 do begin
      grabslit = where(sciobj.slitid EQ slitnum[ii])
      objnum[ii] = sciobj[grabslit[0]].objid
   endfor
endif 

openw, 12, outfil
for ii=0, nslits-1 do begin

   printf,12, slitnum[ii], name[ii], mag[ii], objnum[ii],$
          format='(i3,2x,a14,2x,f5.2,4x,i3,2x)'
endfor   
close, 12

end
