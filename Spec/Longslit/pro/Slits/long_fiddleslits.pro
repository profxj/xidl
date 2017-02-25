;+
; NAME:
;   long_fiddleslits
;
; PURPOSE:
;   Find the location of objects within each slit mask
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;  KEEP -- Slit numbers to keep (starts counting at 1)
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   
; PROCEDURES CALLED:
;   traceset2xy
;   
; REVISION HISTORY:
;   29-Mar-2010  Written by JXP
;-  
;------------------------------------------------------------------------------
pro long_fiddleslits, slit_fil, REMOVE=remove, KEEP=keep, WAVE_FIL=wave_fil
                      

   if N_PARAMS() LT 1 then begin
     print, 'long_fiddleslits, slit_fil [v1.0]' 
     return
   endif

   ;; Read in 
   img = xmrdfits(slit_fil,0,head=hdr)
   tset = xmrdfits(slit_fil,1)
   sz = size(tset.coeff, /dimens)
   sz2 = size(tset.xcorr_coeff, /dimens)
   nslit = sz[1]
   msk = replicate(1B, nslit)

   if keyword_set(KEEP) and keyword_set(REMOVE) then begin
       print, 'long_fiddleslits:  Cannot KEEP and REMOVE'
       return
   endif

   ;; REMOVE
   if keyword_set(REMOVE) then begin
       msk[REMOVE-1] = 0B
       idl_keep = where(msk)
   endif else begin
       if not keyword_set(KEEP) then stop
       idl_keep = keep-1
   endelse
   nkeep = n_elements(idl_keep)

   tset_proto = $
     { func    :    tset[0].FUNC, $
       xmin    :    tset[0].XMIN, $
       xmax    :    tset[0].XMAX, $
       coeff   :    dblarr(3, nkeep),  $
       dims    :    tset[0].DIMS, $
       xcorr_coeff : dblarr(3, sz2[1])  $
     }
   tset_slits = replicate(tset_proto, 2)
   tset_slits[0].COEFF[*, *] = tset[0].COEFF[*, idl_keep]
   tset_slits[1].COEFF[*, *] = tset[1].COEFF[*, idl_keep]
   tset_slits[0].XCORR_COEFF = tset[0].xcorr_coeff
   tset_slits[1].XCORR_COEFF = tset[1].xcorr_coeff

   ;; Make the new slitmask
   slitmask = long_slits2mask(tset_slits, nslit = nslit)

   ;; Write
   splog, 'long_fiddleslits: Writing output file '+slit_fil
   mwrfits, slitmask, slit_fil, hdr, /create
   zeropix = where(slitmask LT 1., nzero)
   mwrfits, tset_slits, slit_fil

   if keyword_set(WAVE_FIL) then begin
       waveimg = xmrdfits(wave_fil,0,head=whdr)
       if nzero NE 0 then waveimg[zeropix] = 0.
       wave2dstr = xmrdfits(wave_fil,1)
       wave1dstr = xmrdfits(wave_fil,2)
       wave2dstr = wave2dstr[idl_keep]
       wave1dstr = wave1dstr[idl_keep]
       ;; Write
       mwrfits, waveimg, wave_fil, hdr, /create 
       mwrfits, wave2dstr, wave_fil 
       mwrfits, wave1dstr, wave_fil 
   endif else begin
       print, 'long_fiddleslits:  You should remake the wavelength image now'
       print, 'long_fiddleslits:  Or you could have had this code do it by giving it the wave image FITS file.  Too late now!' 
   endelse
   
   return
end
;------------------------------------------------------------------------------
