;+ 
; NAME:
; hires_mkascii   
;     Version 1.0
;
; PURPOSE:
;   Convert the binary FITS tables into ASCII files
;
; CALLING SEQUENCE:
;  hires_mkascii, filnm
;
; INPUTS:
;   filnm --  Filename of combined 2D spectrum (output from
;            hires_combspec)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  OUTDIR=  -- Output directory [default:  IFIL/]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Sep-2005 Written by JXP
;-
;------------------------------------------------------------------------------

pro hires_mkascii, filnm, OUTDIR=outdir, EXTRACT=EXTRACT


  ; 
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'hires_mkascii, filnm, OUTDIR= (v1.0)'
      return
  endif 

  if not keyword_set( OUTDIR ) then outdir = 'IFIL/'
;


  if keyword_set(EXTRACT) then begin
      obj = xmrdfits(filnm, 1,/silent)
      nordr = n_elements(obj)

      ;; Exposure number
      ipos = strpos(filnm, 'Obj_')
      if ipos LT 0 then stop else cexp = strmid(filnm,ipos+4,4)

      ;; Loop
      for jj=0L,nordr-1 do begin
          if obj[jj].order LT 100 then $
            co = '0'+strtrim(obj[jj].order,2) $
          else co = strtrim(obj[jj].order,2)
          ;; 
          gd = where(obj[jj].wave GT 0. and obj[jj].var GT 0., ngd)
          if ngd EQ 0 then continue
          id1 = gd[0]
          id2 = gd[ngd-1]
          ;; 
          outfil = outdir+strtrim(obj[jj].field,2)+'_'+cexp+'_'+co+'.asc'
          writecol, outfil, obj[jj].wave[id1:id2], $
                    obj[jj].fx[id1:id2], $
                    sqrt(obj[jj].var[id1:id2]), $
                    FMT='(f10.4,1x,f14.4,1x,f14.4)'
      endfor
  endif else begin
      ;; FSpec format
      obj = xmrdfits(filnm, 1,/silent)
      use = where(obj.phys_ordr NE 0,nordr)
      
      ;; Loop
      for jj=0L,nordr-1 do begin
          ii = use[jj]
          if obj.phys_ordr[ii] LT 100 then $
            co = '0'+strtrim(obj.phys_ordr[ii],2) $
          else co = strtrim(obj.phys_ordr[ii],2)
          ;; 
          gd = where(obj.wave[*,ii] GT 0. and obj.var[*,ii] GT 0., ngd)
          if ngd EQ 0 then continue
          id1 = gd[0]
          id2 = gd[ngd-1]
          ;; 
          outfil = outdir+strtrim(obj.field,2)+'_'+co+'.asc'
          writecol, outfil, obj.wave[id1:id2,ii], $
                    obj.fx[id1:id2,ii], $
                    obj.var[id1:id2,ii], $
                    FMT='(f10.4,1x,f14.4,1x,f14.4)'
      endfor
  endelse
      
  return
end
  
      
