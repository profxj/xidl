;+ 
; NAME:
; hires_writefits   
;     Version 1.0
;
; PURPOSE:
;    Writes one of the standard output files to an ASCII file or to
;    simple FITS output files.
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_writefits, hires, indx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Sep-2005 Written by JXP
;-
;------------------------------------------------------------------------------

pro hires_writefits, hires, indx, SILENT=silent, OUTDIR=outdir, $
                     ASCII=ascii, FITS=fits, PHEAD=phead, AFLUX=aflux, $
                     FSPEC=fspec


  ; 
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_writefits, hires, indx, /silent, /ASCII, /PHEAD, /AFLUX (v1.0)'
      return
  endif 

  if not keyword_set(OBJ_NM) then obj_nm = 'a'
  if not keyword_set( OUTDIR ) then outdir = 'IFIL/'
;

  ;; Loop
  nindx =n_elements(indx)

  if not keyword_set(FSPEC) then begin
      for qq=0L,nindx-1 do begin
          ;;  
          idx = indx[qq]
          objfil = hires_getfil('obj_fil', /name, $
                                FRAME=hires[indx[qq]].frame,$
                                CHIP=hires[indx[qq]].chip, chkfil=chkf)
          if chkf EQ 0 then continue
          ;; Write
          obj = xmrdfits(objfil, 1, /silent)
          
          ;; Header
          finfil = hires_getfil('fin_fil', /name, $
                                FRAME=hires[indx[qq]].frame,$
                                CHIP=hires[indx[qq]].chip, chkfil=chkf)
          head = xheadfits(finfil)
          sxaddpar, head, 'CRPIX1', 1
          sxaddpar, head, 'CTYPE1', 'LINEAR'
          sxaddpar, head, 'DC-FLAG', 1
          sxaddpar, head, 'BITPIX', -32
          sxaddpar, head, 'NAXIS', 1
          sxdelpar, head, 'NAXIS2'
          ;;
          nordr = n_elements(obj)
          case hires[idx].chip of
              1: cc = 'B'
              2: cc = 'G'
              3: cc = 'R'
              else: stop
          endcase
          
          dum = hires[idx].frame
          cframe = ''
          while( DUM LT 1000 ) do begin
              dum = dum*10
              cframe = cframe+'0'
          endwhile
          cframe = cframe+strtrim(hires[idx].frame,2)
          
          for ii=0L,nordr-1 do begin
              if obj[ii].order LT 100 then co = '0'+strtrim(obj[ii].order,2) $
              else co = strtrim(obj[ii].order,2)
              if obj[ii].npix LE 1 then continue
              sxaddpar, head, 'NAXIS1', n_elements(obj[ii].npix)
              if keyword_set( FITS) then begin
                  ;; Wavelengths
                  outwv = outdir+hires[idx].obj+'_'+cframe+'_'+cc+co+'_wv.fits'
                  mwrfits, obj[ii].wave[0:obj[ii].npix-1], outwv, head, /create
                  
                  ;; Flux
                  outfx = outdir+hires[idx].obj+'_'+cframe+'_'+cc+co+'_fx.fits'
                  mwrfits, obj[ii].fx[0:obj[ii].npix-1], outfx, head, /create
                  
                  ;; Sigma
                  outer = outdir+hires[idx].obj+'_'+cframe+'_'+cc+co+'_er.fits'
                  mwrfits, obj[ii].sig[0:obj[ii].npix-1], outer, head, /create
              endif
              if keyword_set( ASCII ) then begin

                  outfil = outdir+hires[idx].obj+'_'+cframe+'_'+cc+co+'.asc'
                  writecol, outfil, obj[ii].wave[0:obj[ii].npix-1], $
                            obj[ii].fx[0:obj[ii].npix-1], $
                            obj[ii].var[0:obj[ii].npix-1], $
                            FMT='(f10.4,1x,f14.4,1x,f10.4)'
              endif
              if keyword_set( AFLUX ) then begin
                  outfil = outdir+hires[idx].obj+'_'+cframe+'_'+cc+co+'.afx'
                  writecol, outfil, obj[ii].wave[0:obj[ii].npix-1], $
                            obj[ii].flux[0:obj[ii].npix-1], $
                            obj[ii].sig[0:obj[ii].npix-1], $
                            FMT='(f10.4,1x,f14.4,1x,f10.4)'
              endif
          endfor
          if keyword_set( PHEAD ) then begin
              outfil = outdir+hires[idx].obj+'_'+cframe+'_'+cc+'.hdr'
              writecol, outfil, head, replicate(' ', n_elements(head))
          endif
      endfor
  endif else begin
      for qq=0L,nindx-1 do begin
          ;;  
          idx = indx[qq]
          subfil = strcompress(strtrim(hires[idx].Obj,2),/remove_all)+obj_nm
          case hires[idx].chip of 
              1: clrc = '_B' 
              2: clrc = '_G'
              3: clrc = '_R'
          endcase
          objfil = 'FSpec/'+subfil+clrc+'.fits'

          ;; Read
          obj = xmrdfits(objfil, 1, /silent)
          
          ;; Header
          finfil = hires_getfil('fin_fil', /name, $
                                FRAME=hires[indx[qq]].frame,$
                                CHIP=hires[indx[qq]].chip, chkfil=chkf)
          head = xheadfits(finfil)
          sxaddpar, head, 'CRPIX1', 1
          sxaddpar, head, 'CTYPE1', 'LINEAR'
          sxaddpar, head, 'DC-FLAG', 1
          sxaddpar, head, 'BITPIX', -32
          sxaddpar, head, 'NAXIS', 1
          sxdelpar, head, 'NAXIS2'
          ;;
          gdo = where(obj.phys_ordr GT 0, nordr)

          case hires[idx].chip of
              1: cc = 'B'
              2: cc = 'G'
              3: cc = 'R'
              else: stop
          endcase
          
          for jj=0L,nordr-1 do begin
              ii = gdo[jj]
              if obj.phys_ordr[ii] LT 100 then $
                co = '0'+strtrim(obj.phys_ordr[ii],2) $
              else co = strtrim(obj.phys_ordr,2)

              gdp = where(obj.wave[*,ii] GT 0., npix)
              if npix LE 1 then continue
              sxaddpar, head, 'NAXIS1', npix
              if keyword_set( FITS) then begin
                  ;; Wavelengths
                  crval1 = alog10(obj.wave[gdp[0],ii])
                  velpix = 1.3d * hires[idx].rowbin
                  cdelt1 = alog10(1.0d + velpix / 299792.458d)
                  sxaddpar, head, 'CRVAL1', crval1
                  sxaddpar, head, 'CDELT1', cdelt1
                  
                  ;; Flux
                  outfx = outdir+hires[idx].obj+'_'+cc+co+'_fx.fits'
                  mwrfits, obj.fx[gdp,ii], outfx, head, /create
                  
                  ;; Var
                  outer = outdir+hires[idx].obj+'_'+cc+co+'_var.fits'
                  mwrfits, obj.var[gdp,ii], outer, head, /create
              endif
              if keyword_set( ASCII ) then begin
                  IF NOT KEYWORD_SET(cframe) THEN $
                    outfil = outdir+hires[idx].obj+'_'+cc+co+'.asc' $
                  ELSE outfil = outdir+hires[idx].obj+'_'+cframe+'_'+cc+co+'.asc'
                  writecol, outfil, obj[ii].wave[0:obj[ii].npix-1], $
                            obj[ii].fx[0:obj[ii].npix-1], $
                            obj[ii].var[0:obj[ii].npix-1], $
                            FMT='(f10.4,1x,f14.4,1x,f10.4)'
              endif
              if keyword_set( AFLUX ) then begin
                  outfil = outdir+hires[idx].obj+'_'+cframe+'_'+cc+co+'.afx'
                  writecol, outfil, obj[ii].wave[0:obj[ii].npix-1], $
                            obj[ii].flux[0:obj[ii].npix-1], $
                            obj[ii].sig[0:obj[ii].npix-1], $
                            FMT='(f10.4,1x,f14.4,1x,f10.4)'
              endif
          endfor
          if keyword_set( PHEAD ) then begin
              outfil = outdir+hires[idx].obj+'_'+cframe+'_'+cc+'.hdr'
              writecol, outfil, head, replicate(' ', n_elements(head))
          endif
      endfor
  endelse
      
  return
end
  
      
