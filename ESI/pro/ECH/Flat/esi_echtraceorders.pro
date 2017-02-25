
;+ 
; NAME:
; esi_echtraceorders
;     Version 1.1
;
; PURPOSE:
;    Traces the order edges
;
; CALLING SEQUENCE:
;   
;  esi_echtraceorders, esi, /CHK,/IFLAT
;
; INPUTS:
;   esi     -  ESI structure
;   slit    -  Slit width  (0.5, 0.75, 1.0)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /CHK      - Display final image and a few other steps
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Only set for 1x1 binning -- Many 'hard wired' numbers
;
; EXAMPLES:
;   esi_echfltsct, esi, 0.75, /CHK
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   11-Aug-2002 Written by JXP
;   01-Feb-2003 Polished (JXP)
;   17-Feb-2004 Added kludge for order #2
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echtraceorders, esi, slit, CHK = chk, IFLAT = IFLAT, SEDG_FIL=sedg_fil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echtraceorders, esi, slit, /IFLAT, /CHK, [v2.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( CBIN ) then cbin = 1
  if not keyword_set( RBIN ) then rbin = 1
  if keyword_set(IFLAT) then ftype = 'I' else ftype = 'D'
    
  ;; Open Flat
  if ftype EQ 'I' then mode = 0 else mode = 1
  flatfil = esi_getfil('echflat_fil', mode, SLIT = slit, $
                       cbin = cbin, rbin = rbin, /name)
  print, 'esi_echtraceorders: Opening flat ', flatfil
  if x_chkfil(flatfil+'*') EQ 0 then begin
      print, 'esi_echtraceorders: File ', flatfil $
             , ' does not exist.  Returning!'
      return
  endif
  flat = xmrdfits(flatfil, 0, head, /silent)
  sz_flat = size(flat, /dimensions)
  nx = sz_flat[0]
  ny = sz_flat[1] 
  
  slit_edg = fltarr(sz_flat[1], 10, 2)
  tset_slits = esi_traceorders(flat)
  traceset2xy, tset_slits[0], rows, left_edge
  traceset2xy, tset_slits[1], rows, right_edge
  slit_edg[*, *, 0] = left_edge
  slit_edg[*, *, 1] = right_edge
  
  ;; CHK (Plot slit edges)
  if keyword_set( CHK ) then begin
      tmp = flat
      for i=0,9 do begin
          for j=0,1 do begin
              rnd_trc2 = round(slit_edg[*,i,j])
              trc_msk = rnd_trc2 + lindgen(sz_flat[1])*sz_flat[0]
              tmp[trc_msk] = -10000
              if keyword_set(IFU) then begin
                  if j EQ 0 then begin
                      rnd_trc2 = round(slit_edg[*,i,0]+ $
                                       (slit_edg[*,i,1]-slit_edg[*,i,0])*0.53)
                      trc_msk = rnd_trc2 + lindgen(sz_flat[1])*sz_flat[0]
                      tmp[trc_msk] = -10000
                      rnd_trc2 = round(slit_edg[*,i,0]+ $
                                       (slit_edg[*,i,1]-slit_edg[*,i,0])*0.57)
                      trc_msk = rnd_trc2 + lindgen(sz_flat[1])*sz_flat[0]
                      tmp[trc_msk] = -10000
                  endif
              endif
          endfor
      endfor
      xatv, tmp, /block, min = -10.0, max = 4.0d4
  endif
  ;; Edge file
  sedg_fil = esi_getfil('sedg_fil', SLIT=slit, cbin=cbin, rbin=rbin, /name)
  print, 'esi_echtraceorders: Writing slit edge info to: ', sedg_fil
  mwrfits, slit_edg, sedg_fil, /create
  mwrfits, tset_slits, sedg_fil
  print, 'esi_echtraceorders: All done!'

  return 
end



