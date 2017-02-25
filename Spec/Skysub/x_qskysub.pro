;+ 
; NAME:
; x_qskysub   
;    Version 1.1
;
; PURPOSE:
;  Sky subtract column by column using a POLY fit.
;
; CALLING SEQUENCE:
;  x_qskysub, data_fil, slit_fil, slit, subimg, pix
;
; INPUTS:
;   data_fil    - Name of Image file
;   slit_fil    - File for slit structure
;   slit        - Slit to sky subtract
;
; RETURNS:
;   subimg  -- Sub area of the img array corresponding to the slit
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
; NWIMG -- Sky subtracted image
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_qskysub, data_fil, slit_fil, slit, subimg, pix, NWIMG=nwimg


;  Error catching
  if  N_params() LT 5  then begin 
    print,'Syntax - ' + $
             'x_qskysub, data_fil, slit_fil, slit, subimg, pix [v1.1]'
    return
  endif 


;  Optional Keywords

;  Read in the data and wave info
   flux = mrdfits(data_fil)
   wave = mrdfits(data_fil,2)
   sz = size(flux, /dimensions)

;  Read in slitstruct
   slitstr = mrdfits(slit_fil, 1, STRUCTYP='mslitstrct', /silent)

;  Find all pixels in that slit
   msk = bytarr(sz[0],sz[1])
   for i=0L, sz[0]-1 do begin
       yedg = round(slitstr[slit].yedg_orig[i,*])
       msk[i,yedg[0]+1:yedg[1]-1] = 1
   endfor

;  Create fit structure
   fitstr = { fitstrct }
   fitstr.func = 'POLY'
   fitstr.nord = 2
   fitstr.hsig = 3.
   fitstr.lsig = 5.
   fitstr.maxrej = round(0.1*(slitstr[slit].yedg_flt[1]-$
                              slitstr[slit].yedg_flt[0]))
   fitstr.niter = 3
   fitstr.minpt = 3

;  Loop and Fit
   nwimg = fltarr(sz[0], sz[1])
   for i=0L,sz[0]-1 do begin
       gdpix = where(msk[i,*] EQ 1)
       sky = x_fitrej(wave[i,gdpix], flux[i,gdpix], FITSTR=fitstr)
       nwimg[i,gdpix] = flux[i,gdpix] - sky
   endfor

;  Create outputs   
   pix = where(msk EQ 1)
   subimg = nwimg[pix]

   return
end

   
   
