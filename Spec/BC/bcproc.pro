;+ 
; NAME:
; bcproc   
;   Version 1.0
;
; PURPOSE:
;    Processes a BC Spectrograph image:  Subtract OV, BIAS, divide
;    flat
;
; CALLING SEQUENCE:
;   
;   bcproc, filename, image, [invvar, hdr,] bias=, flat=
;
; INPUTS:
;   filename   - String
;
; RETURNS:
;   image      - Processed image
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   bias       - Bias image
;   flat       - Flat image
;
; OPTIONAL OUTPUTS:
;   invvar     - Inverse variance
;   hdr        - Header of the input image
;
; COMMENTS:
;
; EXAMPLES:
;   bcproc, science, image, invvar, hdr, bias=bias, flat=flat
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Aug-2001 Written by SMB
;   10-Oct-2001 Minor modifications JXP
;-
;------------------------------------------------------------------------------

pro bcproc, filename, image, invvar, hdr, bias=bias, flat=flat

;
    if  N_params() LT 4  then begin 
        print,'Syntax - ' + $
          'bcproc, file, image, [invvar, hdr,] bias=, flat= (v1.0)'
        return
    endif 

    tt = mrdfits(filename,0,hdr)

    if (size(tt))[0] NE 2 then return

    ncol = (size(tt))[1]
    nrow = (size(tt))[2]

    readoutvariance = 0.0

    if ncol GT 1026 then begin
      nover = lindgen(ncol - 1024) + 1024L
      ncol = 1022L
      djs_iterstat, tt[nover,*], sigma=readoutnoise
      readoutvariance = readoutnoise^2
      overs = replicate(1,ncol) # djs_median(tt[nover,*],1) 
      image = tt[0:ncol-1,*] - overs
    endif else image = tt

    invvar = 1.0/(abs(image) + (image EQ 0) + readoutvariance)

    if keyword_set(bias) then $
      if n_elements(bias) EQ n_elements(image) EQ 0 then $
         image = image - bias

    if keyword_set(flat) then $
       divideflat, image, invvar, flat, minval=0.3

    return
end
