;+ 
; NAME:
; x_qskysub   
;    Version 1.0
;
; PURPOSE:
;    Given the slitstr and the map, find slit positions in the
;    original image
;
; CALLING SEQUENCE:
;   
;   x_qskysub, slitstr, map
;
; INPUTS:
;   slitstr     - Slit structure
;   map         - y-distortion map (fits is ok)
;
; RETURNS:
;
; OUTPUTS:
;   Updates slitstr for original positions
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_qskysub, slitstr, map
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_qskysub, data_fil, slit_fil, slit, subimg, pix


;  Error catching
  if  N_params() LT 5  then begin 
    print,'Syntax - ' + $
             'x_qskysub, data_fil, slit_fil, slit, subimg, pix [v1.0]'
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

   
   
