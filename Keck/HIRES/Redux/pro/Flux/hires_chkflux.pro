;+ 
; NAME:
; hires_chkflux
;     Version 1.1
;
; PURPOSE:
;  Creates a series of plots to show the end result of fluxing.  In
;  particular one can check the fluxing in the order overlap regions.
;
; CALLING SEQUENCE:
;  hires_chkflux, hires, setup, obj_id, chip, [exp]
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
;   FIL -  File name of standard star for fluxing (default:
;               $MIKE_DIR/pro/Std/Archive/sens_blue#.fits)
;   /COMB -
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_chkflux, hires, setup, obj_id, chip, /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-Sep-2005 JXP
;-
;------------------------------------------------------------------------------

pro hires_chkflux, hires, setup, obj_id, chip, exp, PRINT=print, FIL=fil, $
                   COMB=comb

  if  N_params() LT 4 and not keyword_set(FIL) then begin 
      print,'Syntax - ' + $
        'hires_chkflux, hires, setup, obj_id, chip, [exp], /PRINT [v1.1]'
      return
  endif

  ;; Optional Keywords
  if not keyword_set(CSZ) then csz = 3.

  if not keyword_set(PRINT) then begin
      device, get_screen_size=ssz
      if not keyword_set( XSIZE ) then    xsize = ssz[0]-200
      if not keyword_set( YSIZE ) then    ysize = ssz[1]-200
  endif else csz = 2.
  if not keyword_set(NPY) then npy = 3
  if not keyword_set(NSPEC) then nspec = 4


; Grab exposure
  if not keyword_set( FIL ) then begin
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
      allexp = 0
      nexp = 1
  endelse


  ;; Loop on exposures
  for qq=0L,nexp-1 do begin
      idx = allexp[qq]
      ;; Print?
      if keyword_set(PRINT) then begin
          if not keyword_set(FIL) then $
            psfile = hires_getfil('qa_chkflux',setup, frame=hires[idx].frame, $
                                  chip=chip,/name) $
          else psfile = fil+'.ps'
          x_psopen, PSFILE, /maxs
      endif
      !p.multi=[0,1,npy]
      clr = getcolor(/load)

      
      if not keyword_set(PRINT) then begin
          if keyword_set(FIL) then obj = fil else obj = hires[idx].Obj
          window, 1, title=Obj+'  Exposure = '+strtrim(qq), $
                  XSIZE=xsize, YSIZE=ysize
      endif

      ;; Loop on orders
      if not keyword_set(FIL) then objfil = hires[idx].obj_fil else objfil = fil
      objstr = xmrdfits(objfil, 1, /silent)

      if keyword_set(COMB) then begin
          sz = size(objstr.fx,/dimensions)
          sumv = total(objstr.var, 1)
          gdo = where(objstr.phys_ordr GT 0 AND $
                      sumv GT 0., nord)
          tmp = { $
                   wave: dblarr(sz[0]), $
                   flux: fltarr(sz[0]), $
                   sig: fltarr(sz[0]), $
                   order: 0 $
          }
          newstr = replicate(tmp, nord)
          newstr.wave = objstr.wave[*,gdo]
          newstr.flux = objstr.fx[*,gdo]
          newstr.sig = sqrt(objstr.var[*,gdo])
          newstr.order = objstr.phys_ordr[gdo]
          objstr = newstr
      endif

      nord = n_elements(objstr)
      nplt = nord / (npy*nspec) + ((nord mod (npy*nspec)) NE 0)

      ;;  Plot
      ntot = 0L
      for nn=0L,nplt-1 do begin
          for ii=0L,npy-1 do begin
              ;; Indices
              obji = ntot + lindgen(nspec)
              obji = obji < (nord-1)
              gd = where(objstr[obji].sig GT 0.)
              ;; Ranges
              mdf = median((objstr[obji].flux)[gd])
              mnwv = min((objstr[obji].wave)[gd],max=mxwv)
              plot, [0], [0],  xrange=[mnwv,mxwv], background=clr.white, $ 
                    color=clr.black,$ ; xtickn=xspaces, $
                    xtitle='Wavelength', xstyle=1, ystyle=1, $
                    charsize=csz, /nodata, $
                    xmargin=xrmg, ymargin=[0,0.], yrange=[-0.2*mdf,mdf*3]

              ;; Plot
              for jj=0L,nspec-1 do begin
                  case jj of
                      0: cclr = clr.red
                      1: cclr = clr.blue
                      2: cclr = clr.black
                      3: cclr = clr.green
                      4: cclr = clr.yellow
                      else: stop
                  endcase

                  gdw = where(objstr[obji[jj]].wave GT 0.)
                  oplot, objstr[obji[jj]].wave[gdw], thick=2, $
                         objstr[obji[jj]].flux[gdw], color=cclr, psym=10
              endfor
              ntot = ntot + nspec
          endfor
          if not keyword_set(PRINT) then stop 
      endfor
      if keyword_set( PRINT ) then begin
          x_psclose 
          spawn, 'gzip -f '+psfile
      endif
  endfor

  !p.multi=0
  print, 'hires_chkflux: All done'
  return

end

