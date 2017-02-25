;+
; 
; NAME:
; xdimg_badpix
;    Version 1.1
;
; PURPOSE:
;    Creates a bad pixel mask given a structure which contains either a flat or dark image.
;
; CALLING SEQUENCE:
;   xdimg_badpix, struct
;
; INPUTS:
;   struct -- dimg_strct defining the images of interest
;
; RETURNS:
;
; OUTPUTS:
;   OUTFIL - Name of bad pixel mask file [default is e.g. 'PixMask/BadPixMask_w.fits']
;            NOTE: the "w" indicates NIRC2 wide camera and the "n" indicates narrow camera.
;
; OPTIONAL KEYWORDS:
;   OUTFIL - Name of bad pixel mask to be created, including directory.
;   CCD - Name of the ccd for which a bad pixel mask should be made.  ('NIRC2W' or 'NIRC2N')
;   DEADPERCENT - percentage of saturation to be identified as a dead pixel in a well-exposed flat frame. (default is 10%)
;   HOTPERCENT - percentage of saturation to be identified as a hot pixel in a dark frame. (default is 50%)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xdimg_badpix, nght1_strct
;
; PROCEDURES/FUNCTIONS CALLED:
;  XDIMG_DELCODIV
;
; REVISION HISTORY:
;   11-June-2007 Written by LKP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_badpix, struct, OUTFIL=outfil, CCD=ccd, DEADPERCENT=deadpercent, HOTPERCENT=hotpercent
  
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'xdimg_badpix, struct, OUTFIL=outfil, CCD=ccd, DEADPERCENT=deadpercent, HOTPERCENT=hotpercent (v1.1)'
      return
  endif 

  if not keyword_set( DEADPERCENT ) then deadpercent = 0.1
  if not keyword_set( HOTPERCENT ) then hotpercent = 0.5

  use = where(struct.flg_anly NE 0, usenum)
  xdim = struct[use[0]].naxis1
  ydim = struct[use[0]].naxis2
  mask = fltarr(xdim, ydim)
  mask[*,*]=1.


;  Find the domeflats
  if keyword_set(CCD) then dflts = where(struct.type EQ 'DFT' AND struct.flg_anly NE 0 AND struct.ccd EQ ccd, ndflt) $
                      else dflts = where(struct.type EQ 'DFT' AND struct.flg_anly NE 0, ndflt)
  if ndflt eq 0 then begin
     print, 'Not using domeflats to make a badpixel mask.'  
  endif else begin

     ; Loop over all domeflats available.
     for i=0,ndflt-1 do begin

        ; Preprocessing -- Divide flat by coadds to compare counts correctly.
        ; Don't worry about dark subtracting here.
        flatfil = strtrim(struct[dflts[i]].rootpth,2) + strtrim(struct[dflts[i]].img_root,2)
        flat = readfits(flatfil, /silent)
        flat = flat/float(struct[dflts[i]].coadds)

        ; Find dead pixels in a well-exposed flat frame.
        med = median(flat)
        if med lt 4000 then begin
           ndead = 0
           dead = 0
        endif else begin
           dead = where(flat lt struct[dflts[i]].satur * deadpercent, ndead)
        endelse

        ; Update the mask with the dead pixels in each of the flats.
        if ndead gt 0 then mask[dead]=0.

     endfor
     
  endelse
  

;  Find the twiflats
  if keyword_set(CCD) then twiflts = where(struct.type EQ 'TWI' AND struct.flg_anly NE 0 AND struct.ccd EQ ccd, ntwiflt) $
                      else twiflts = where(struct.type EQ 'TWI' AND struct.flg_anly NE 0, ntwiflt)
  if ntwiflt eq 0 then begin
     print, 'Not using twilight flats to make a badpixel mask.'  
  endif else begin

     ; Loop over all twilight flats available.
     for i=0,ntwiflt-1 do begin

        ; Preprocessing -- Divide flat by coadds.
        ; Don't worry about dark subtracting here.
        flatfil = strtrim(struct[twiflts[i]].rootpth,2) + strtrim(struct[twiflts[i]].img_root,2)
        flat = readfits(flatfil, /silent)
        flat = flat/float(struct[twiflts[i]].coadds)
        
        ; Find dead pixels in a well-exposed twilight frame.
        med = median(flat)
        if med lt 4000 then begin
           ndead = 0
           dead = 0
        endif else begin
           dead = where(flat lt struct[dflts[i]].satur * deadpercent, ndead)
        endelse

        ; Update the mask with the dead pixels in each of the flats.
        if ndead gt 0 then mask[dead]=0.

     endfor
     
  endelse

  
  ; Find dark frames and look for hot pixels.
  if keyword_set(CCD) then darkims = where(struct.type EQ 'DRK' AND struct.flg_anly NE 0 AND struct.ccd EQ ccd, ndark) $
                      else darkims = where(struct.type EQ 'DRK' AND struct.flg_anly NE 0, ndark)
  if ndark eq 0 then begin
     print, 'Not using dark frame to make a badpixel mask.'
  endif else begin

     ; Loop over all dark frames available.
     for i=0, ndark-1 do begin
        
        ; Preprocessing -- Divide by coadds
        darkfil = strtrim(struct[darkims[i]].rootpth,2) + strtrim(struct[darkims[i]].img_root,2)
        dark = readfits(darkfil, /silent)
        dark = dark/float(struct[darkims[i]].coadds)
        
        ; Find hot pixels in a dark frame.
        hot = where(dark gt struct[darkims[i]].satur * hotpercent, nhot)
        
        ; Update the mask with the hot pixels in each of the darks.
        if nhot gt 0 then mask[hot]=0.

     endfor
  endelse
  

if ndflt EQ 0 and ndark EQ 0 and ntwiflt EQ 0 then begin
   print, 'Cannot make bad pixel mask.  Have neither darks nor flats to work with.'
   return
endif


; Save the mask
  if keyword_set(CCD) then begin
     CASE ccd of
        'NIRC2W': begin
           if not keyword_set( OUTFIL ) then outfil = 'PixMask/BadPixMask_w.fits'
        end
        'NIRC2N': begin
           if not keyword_set( OUTFIL ) then outfil = 'PixMask/BadPixMask_n.fits'
        end
     end
  endif else begin
     if not keyword_set( OUTFIL ) then outfil = 'PixMask/BadPixMask.fits'
  endelse
  writefits, strtrim(outfil,2), mask


; Resave the updated structure so that you can pick up where you left off easily.
  mwrfits, struct, 'struct.fits', /create

  print, 'All done with bad pixel mask!'

end
