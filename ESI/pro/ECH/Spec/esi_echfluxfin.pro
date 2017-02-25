;+ 
; NAME:
; esi_echfluxfin  
;     Version 1.1
;
; PURPOSE:
;    Fluxes the spectrum using a standard star calibration.  Default
;    is to use the one in CALIBS which is probably good enough for
;    relative fluxing.  Definitely not good enough for absolute.
;
; CALLING SEQUENCE:
;   
;  esi_echfluxfin, esi, obj_id, [fluxfil], /CLOBBER, OBJ_NM=, /STD
;
; INPUTS:
;   esi       -  ESI structure
;   obj_id    -  Object ID  (e.g. 0L, 1L, etc)
;   [fluxfil] -  File name of standard star for fluxing (default:
;               $XIDL_DIR/ESI/CALIBS/ECH_FLUX_050.fits)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;    ORDRS=    - Orders to flux (default: [0L,9L])
;   SPECFIL=  - Useful for fluxing files which resulted from a
;               combination of multiple nights
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echfluxfin, esi, obj_id, /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Oct-2002 Written by JXP
;   04-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------

pro esi_echfluxfin, esi, obj_id, fluxfil = fluxfil, CLOBBER = clobber, $
                    OBJ_NM = obj_nm, $
                    STD=std, ORDRS=ordrs, SPECFIL=specfil, OLD=old

  if  N_params() LT 2  and NOT keyword_set( SPECFIL ) then begin 
      print,'Syntax - ' + $
        'esi_echfluxfin, esi, obj_id, [fluxfil], /CLOBBER, OBJ_NM=, /STD '
      print, '     ORDRS=, /OLD [v1.1]'
      return
  endif

;  Optional Keywords

  if not keyword_set(ORDRS) then ordrs=[0L,9L]
  if not keyword_set(OBJ_NM) then obj_nm = 'a'
  if not keyword_set( FLUXFIL ) then $
    fluxfil = getenv('ESI_CALIBS')+'/ECH_FLUX_050.fits'

; Grab exposure
  if not keyword_set( SPECFIL ) then begin
      if not keyword_set( STD ) then begin
          allexp = where(esi.type EQ 'OBJ' AND esi.flg_anly NE 0 AND $
                         esi.mode EQ 2 AND esi.obj_id EQ obj_id, nexp)
          if nexp EQ 0 then begin
              print, 'esi_echfluxfin: No objects found!'
              return
          endif
      endif else begin
          allexp = obj_id[0]
          nexp = 1
      endelse

      ;; Specfil
      specfil = 'FSpec/'+strtrim(esi[allexp[0]].Obj,2)+obj_nm+'_ech.fits'
  endif

  ;; Open specfil
  if x_chkfil(specfil+'*') EQ 0 then begin
      print, 'esi_echfluxfin: Spec file doesnt exist! Returning..', specfil
      return
  endif

  x_wrechfspec, spec, specfil, /READ
  if spec.flg_flux EQ 1 AND not keyword_set(CLOBBER) then begin
      print, 'esi_echfluxfin: Already fluxed!'
      return
  endif

  ;; spectra are scaled to zeroth member, so use zeroth exp time
  IF NOT KEYWORD_SET(SPECFIL) THEN time = esi[allexp[0]].EXP $
  ELSE time = spec.TEXP[0]    
  
  ;; Flux func
  if keyword_set( OLD ) then begin
      for qq = ordrs[0], ordrs[1] do begin
          ;; Read bset
          bset = xmrdfits(fluxfil, qq+1, /silent)
          
          ;; Get values
          full_rtio = 10^bspline_valu(spec.wave[*,qq], bset)
          
          spec.fx[*,qq] = spec.fx[*,qq]*full_rtio
          b = where(spec.var[*,qq] GT 0.)
          spec.var[b,qq] = spec.var[b,qq]*(full_rtio[b]^2)
          spec.novar[b, qq] = spec.novar[b, qq]*(full_rtio[b]^2)
          spec.sky[b, qq] = spec.sky[b, qq]*full_rtio
      endfor
  endif else begin
     restore, fluxfil
      for qq=ordrs[0],ordrs[1] do begin
          npix = n_elements(spec.fx[*,qq])
          gpx = lindgen(npix)
          kk = where(ordr_fit EQ (15-qq))
          ;; Apply
          full_rtio = x_calcfit(spec.wave[gpx,qq], $
                                fitstr=tot_fit[kk])
          spec.fx[gpx,qq] = spec.fx[gpx,qq]/full_rtio/time
          b = where(spec.var[gpx,qq] GT 0.)
          spec.var[b, qq] = spec.var[gpx[b], qq] $
            / (full_rtio[b]^2) / (time^2)
          spec.novar[b, qq] = spec.novar[gpx[b], qq] $
            / (full_rtio[b]^2) / (time^2)
          spec.sky[gpx, qq] = spec.sky[gpx, qq]/full_rtio/time
      endfor
  endelse
  spec.flg_flux = 1
  
  print, 'esi_echfluxfin: Writing fluxed file ', specfil
  x_wrechfspec, spec, specfil
  print, 'esi_echfluxfin: All done'

end

