;+ 
; NAME:
; hires_flux
;     Version 1.1
;
; PURPOSE:
;    Fluxes the spectrum using a standard star calibration.  Ideal
;    is to use a sensitivity function made from that night.  It is
;    RARE, however, that the fluxing is better than 10% even in a
;    relative sense.    
;
; CALLING SEQUENCE:
;  hires_flux, hires, setup, obj_id, chip, [exp], FLUXFIL=, /CLOBBER, /STD
;
; INPUTS:
;   setup   -  Setup ID
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
;   chip    -  Blue (1) Green (2) or Red (3)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   FLUXFIL -  File name of sensitivity function
;   /STD    - Extraction should be set for a standard star
;   /BOXCAR - Flux data from the boxcar extraction
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_flux, hires, obj_id, /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-Sep-2005 JXP
;-
;------------------------------------------------------------------------------

pro hires_flux, hires, setup, obj_id, chip, exp, EXTRAP=extrap, $
                FLUXFIL=fluxfil, $
                CLOBBER=clobber, STD=std, BOXCAR=boxcar

  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'hires_flux, hires, setup, obj_id, chip, [exp], FLUXFIL=, ' 
      print, '/CLOBBER, /STD, /BOXCAR, /EXTRAP     [v1.1]'
      return
  endif

;  Optional Keywords


; Grab exposure
  if not keyword_set( STD ) then begin
      allexp = where((hires.type EQ 'OBJ' OR hires.type EQ 'STD') $
                     AND hires.flg_anly NE 0 AND $
                     hires.setup EQ setup AND hires.obj_id EQ obj_id $
                     AND hires.chip EQ chip, nexp)
      if nexp EQ 0 then begin
          print, 'hires_flux: No objects found!'
          return
      endif
      ;;  Exposures
      if n_elements(exp) NE 0 then allexp = allexp[exp]
      nexp = n_elements(allexp)
  endif else begin
      allexp = obj_id[0]
      nexp = 1
  endelse

  if not keyword_set( FLUXFIL ) then begin
      stop  ;; Not recommended to use the defaults right now (JXP: 6/26/05)
      stop  ;; Actually, this is not even an option
      return
;      if chip EQ 1 then begin
;          fluxfil = getenv('MIKE_DIR')+'/pro/Std/Archive/sens_blue'+ $
;            strtrim(hires[allexp[0]].rowbin,2)+'.fits'
;      endif else begin
;          fluxfil = getenv('MIKE_DIR')+'/pro/Std/Archive/sens_red'+ $
;            strtrim(hires[allexp[0]].rowbin,2)+'.fits'
;      endelse
      if x_chkfil(fluxfil) EQ 0 then begin
          print, 'hires_flux:  Archived file', fluxfil, 'does not exist!'
          print, 'hires_flux:  You will need to create your own..'
          return
      endif else print, 'hires_flux:  Using archived solution ', fluxfil
  endif
  restore, fluxfil  

; Flux func
  print, 'hires_flux: Fluxing...'
  for qq=0,nexp-1 do begin
      idx = allexp[qq]
      ;; Open objfil
      objfil = hires[idx].obj_fil
      objstr = xmrdfits(objfil, 1, STRUCTYP='hiresobjstrct', /silent)
;      flux = objstr.fx * 0
;      sig  = objstr.var * 0


      ;; Legendre polynomials
      for jj=0L,n_elements(objstr)-1 do begin
          flg_extrap = 0
          ;; Boxcar
          if keyword_set(BOXCAR) then begin
              objstr[jj].npix = objstr[jj].nrow
              objstr[jj].wave = objstr[jj].box_wv
              objstr[jj].fx = objstr[jj].box_fx
              objstr[jj].var = objstr[jj].box_var
          endif
              
          gpx = lindgen(objstr[jj].npix)

          ;; Match up
          mtch = where(objstr[jj].order EQ ordr_fit,nmt)
          if nmt EQ 0 then begin
              print, 'hires_flux: No match to order = ', objstr[jj].order
              if not keyword_set(EXTRAP) then continue else begin
                  mn = min(abs(ordr_fit-objstr[jj].order),mtch)
                  print, 'hires_flux: Extrapolating... -- Good luck!'
                  flg_extrap = 1
                  owv = where(objstr.order EQ ordr_fit[mtch])
                  ewv = objstr[owv].wave[gpx]
              endelse
          endif else mtch = mtch[0]

          ;; Sensitivity function
;              gds = where(sens_wv[*,mtch] GT 0.)
;              linterp, sens_wv[gds,mtch], sens_func[gds,mtch], $
;                objstr[jj].wave[gpx], full_rtio
          if flg_extrap EQ 0 then $
            full_rtio = x_calcfit(objstr[jj].wave[gpx], fitstr=tot_fit[mtch])$
          else $
            full_rtio = x_calcfit(ewv, fitstr=tot_fit[mtch])

          ;; Apply 
          objstr[jj].flux[gpx] = objstr[jj].fx[gpx]/full_rtio/ $
            hires[idx].exp
          
          b = where(objstr[jj].var[gpx] GT 0.,nb)
          if nb NE 0 then begin
              objstr[jj].sig[gpx[b]] = sqrt(objstr[jj].var[gpx[b]]) $
                / full_rtio[b] / hires[idx].exp
          endif
          objstr[jj].flg_flux = 1
      endfor
 
      print, 'hires_flux: Writing fluxed file ', objfil
      mwrfits, objstr, objfil, /create
      spawn, 'gzip -f '+objfil
  endfor
  
  print, 'hires_flux: All done'
  return

end

