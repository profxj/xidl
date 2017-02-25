;+ 
; NAME:
; x_pltechspec
;  V1.1
;
; PURPOSE:
;    Plot the normalized echelle orders in a pretty ps file.
;    system
; CALLING SEQUENCE:
;   
;   x_echspec, strct_fil, instr_list, vel_fil, NTOT=, CSIZE=,
;   LSIZE=, PSFILE=, XTINT=
;
; INPUTS:
;
; RETURNS:
;  strct_fil -- FITS file for the FUSE abs lin structure
;  instr_list -- List of instrument files
;  vel_fil -- Input file for velocity plot
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  NTOT -- Number of plots per page [default: 16]
;  LSIZE -- Label size [default: 1.8]
;  CSIZE -- Numbering character size [default: 1.8]
;  XTINT -- xtick interval
;
; OPTIONAL OUTPUTS:
;  PSFILE -- Postscript filename
;
; COMMENTS:
;
; EXAMPLES:
;   fig_calcewn, struct, fil_instr
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   12-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro x_pltechspec, dat_fil, DWV=dwv, PSFILE=psfile, INFLG=inflg, YR=yr

; lowzovi_prsdat -- Reads in DLA data to a structure

;  if not keyword_set( DAT_FIL ) then $
;    dat_fil = '/u/xavier/Keck/HIRES/RedData/TiII/HD195965/HD195965_f.fits'
  if not keyword_set( LSIZE ) then lsize = 1.6
  if not keyword_set( LTHICK ) then lthick = 3
  if not keyword_set(PSFILE) then psfile = 'x_pltechspec.ps'
  if not keyword_set( CSIZE ) then csize = 1.8
  if not keyword_set(YR) then yr = [-0.1, 1.1]
  if not keyword_set( DWV ) then dwv=100.
  LSIZE = 2.0

; Open vel_fil
  close, /all

  xmrg = [6,0.5]
  ymrg = [3.5,0.5]

  ;; PSFILE
  if keyword_set( PSFILE ) then x_psopen, psfile, /maxs
  !p.multi=[0,1,2,0,1]
  clr = getcolor(/load)

  ;; Read data
  fx = x_readspec(dat_fil, wav=wave, sig=sig, NPIX=npix, INFLG=inflg)
  mnwv = min(wave, max=mxwv)
  nplt = fix((mxwv-mnwv)/dwv) + 1

  for ss=0,nplt-1 do begin 
     xr = mnwv + [dwv*ss, dwv*(ss+1)]

     plot, wave, fx, xrange=xr, yrange=yr, $
           xmargin=xmrg, ymargin=ymrg, /NODATA, $
           charsize=csize, psym=10, background=clr.white, color=clr.black, $
           xstyle=1, ystyle=1, thick=5, ytitle='Normalized Flux', $
           xtitle='Wavelength (Ang)'
  
     ;; Zero line
     oplot, xr, [0., 0.], color=clr.green, linesty=1, thick=1
     oplot, xr, [1., 1.], color=clr.red, linesty=1, thick=2

     ;; Plot
     oplot, wave, fx, color=clr.black, psym=10, thick=3

  endfor
  if keyword_set( PSFILE ) then x_psclose
  !p.multi=[0,1,1]

  return
end
