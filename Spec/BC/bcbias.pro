;+ 
; NAME:
; bcbias   Version 1.1
;
; PURPOSE:
;    Creates a Bias and Flat from a list of BC images
;
; CALLING SEQUENCE:
;   
;   bcbias, files, [biasnm, flatnm]
;
; INPUTS:
;   files      - Vector of strings
;
; RETURNS:
;
; OUTPUTS:
;   Creates a Bias file in 'Bias/'
;   Creates a Flat file in 'Flat/'
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   biasnm     - Name of Bias image
;   flatnm     - Name of Flat
;
; COMMENTS:
;
; EXAMPLES:
;   bcbias, files
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   10-Oct-2001 Minor modifications JXP
;   Aug-2001 Written by SMB
;-
;------------------------------------------------------------------------------

pro bcbias, filenames, biasname, flatname

  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'bcbias, files, [biasname, flatname] (v1.1)'
      return
  endif 

   nfile = n_elements(filenames)
   checkbias = replicate( { ccdpicno : 0L, ncol : 0L, nrow : 0L, $
        utdate : '', utstart: '', object : '', $
        bias : 0L, flat: 0L }, nfile)
   nu = n_elements(u)

   for i=0,nfile - 1 do begin
     hdr = headfits(filenames[i])
     checkbias[i].ccdpicno = sxpar(hdr,'CCDPICNO')
     checkbias[i].ncol = sxpar(hdr,'NAXIS1')
     checkbias[i].nrow = sxpar(hdr,'NAXIS2')
     checkbias[i].utdate = sxpar(hdr,'DATE-OBS')
     checkbias[i].utstart = sxpar(hdr,'UTSTART')
     checkbias[i].object  = sxpar(hdr,'OBJECT')

     checkbias[i].bias = $
         (strpos(strupcase(sxpar(hdr,'OBJECT')),'BIAS') GT -1)  OR $
         (strpos(strupcase(sxpar(hdr,'OBJECT')),'ZERO') GT -1)  OR $
         (strpos(strupcase(sxpar(hdr,'IMAGETYP')),'BIAS') GT -1)  OR $
         (strpos(strupcase(sxpar(hdr,'IMAGETYP')),'ZERO') GT -1) 

     checkbias[i].flat = $
         (strpos(strupcase(sxpar(hdr,'IMAGETYP')),'FLAT') GT -1)  OR $
         (strpos(strupcase(sxpar(hdr,'OBJECT')),'FLAT') GT -1)  
   endfor

   sortstring = checkbias.utdate +'-' +strtrim(string(checkbias.nrow),2)
   s = sort(sortstring)
   u = s[uniq(sortstring[s])]
   nu = n_elements(u)

; Bias First

   a = findfile('Bias/..', count=count)
   if count EQ 0 then file_mkdir, 'Bias'

   for i=0, nu-1 do begin

     here = where(sortstring EQ sortstring[u[i]] AND $
                  checkbias.bias GT 0, nhere)
 
     if nhere GE 3 then begin
       biasname = 'Bias/bcbias-'+sortstring[u[i]]+'.fits'

       bcproc, filenames[here[0]], biashold

       for j=1, nhere-1 do begin
         bcproc, filenames[here[j]], image
         biashold = [[[temporary(biashold)]],[[image]]]
       endfor

       bias = djs_median(biashold,3)
       mwrfits, bias, biasname, /create
     endif
 
   endfor

;  Now do flats:
   a = findfile('Flat/..', count=count)
   if count EQ 0 then file_mkdir, 'Flat'

   u = uniq(checkbias.object)
   nu = n_elements(u)

   for i=0, nu-1 do begin

     here = where(checkbias.object EQ checkbias[u[i]].object AND $
                  checkbias.flat GT 0, nhere)
 
     if nhere GE 3 then begin

       flatname = 'Flat/bcflat-'+fileandpath(filenames[here[0]])

       bcproc, filenames[here[0]], flathold 

       for j=1, nhere-1 do begin
         bcproc, filenames[here[j]], image
         flathold = [[[temporary(flathold)]],[[image]]]
       endfor

       flat = djs_median(flathold,3)
       mwrfits, flat, flatname, /create
     endif
   endfor

return
end
