;+ 
; NAME:
; fuse_velplt
;  V1.1
;
; PURPOSE:
;    Given a list of DLA base files, fill up the structure ;
; CALLING SEQUENCE:
;   
;   lowzovi_prsdat, stucture, filename
;
; INPUTS:
;
; RETURNS:
;   structure      - IDL structure
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  LIST - File
;  ION - Input ionic column densities
;  NOELM - Supress inputting Elemental values
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   fuse_calcewn, struct, fil_instr
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   12-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro fuse_velplt, strct_fil, instr_list, vel_fil, NTOT=NTOT, CSIZE=csize, $
                 LSIZE=lsize, PSFILE=psfile, XTINT=xtint

; lowzovi_prsdat -- Reads in DLA data to a structure

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'fuse_velplt, strct, instr_list, vel_fil (v1.0)' 
    return
  endif 
;  
  if not keyword_set( NTOT ) then ntot = 16L
  if not keyword_set( LSIZE ) then lsize = 1.8
;  if not keyword_set( PSFIL ) then psfil = 'tmp.ps'

;  Read instrument file list
  readcol, instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f'

  nlist = n_elements(instr_fil)
  if nlist LT 7 then stop

;  Open structure
  strct = xmrdfits(strct_fil, 1, /silent)

; Open vel_fil
  close, /all
  openr, 11, vel_fil
  readf, 11, zabs, FORMAT='(f12.7)'
  readf, 11, vmin,vmax
  vmnx = [vmin,vmax]
  nlin = 0
  readf, 11, nlin

; PLOT
  if nlin LE ntot then begin
      if nlin GT 8 then begin
          npx = 2 
          npy = nlin/2 + (nlin MOD 2)
          if not keyword_set( CSIZE ) then csize = 3.0
          xmrg = [2,0]
      endif else begin
          npx = 1
          npy = nlin
          xmrg = [3,0]
          if not keyword_set( CSIZE ) then csize = 3.7
          LSIZE = 2.2
      endelse
  endif else begin
      npx = 2
      npy = ntot/2
      xmrg = [2,0]
      if not keyword_set( CSIZE ) then csize = 3.0
  endelse

  if keyword_set( PSFILE ) then begin
      device, decompose=0
      ps_open, filename=PSFILE, /color, bpp=8, /portrait
      !p.thick = 5
      !p.charthick = 3
  endif

  !p.multi=[0,npx,npy,0,1]
  clr = getcolor(/load)

; Vel strct


  wrest = 0.d
  svinst = 0
  nblnd = 0

; LOOP
  for ii=0L,nlin-1 do begin
      readf, 11, wrest, ymin, ymax, nblnd
      ymnx = [ymin,ymax]

      ;; Get instrument value
      a = where(abs(strct.zabs-zabs) LT 0.001 AND $
                abs(strct.wrest-wrest) LT 0.003, na)
      if na EQ 0 then begin
          print, 'fuse_velplt: No wave! ', wrest
          stop
      endif else a = a[0]
      lgv = alog(double(strct[a].instr)) / alog(2)
      instr = fix(lgv+0.00001)
;      if lgv - jj LT 0.00001 then instr = jj else

      ;; Read data
      if instr NE svinst then begin
          svinst = instr
          ;;  Note that STIS starts at instrument 8
          if instr LE 6 then $
            fx = x_readspec(instr_fil[instr], SIG=sig, wav=wave, $
                            NPIX=npix, inflg=3)$
          else begin  ; STIS
              spos = strpos(instr_fil[instr], 'f.fits')
              sig_fil = strmid(instr_fil[instr], 0, spos)+'e.fits'
              fx = x_readspec(instr_fil[instr], SIG=sig, wav=wave, NPIX=npix, $
                              fil_sig=sig_fil, inflg=0)
          endelse
          ;; Sort?
          srt = sort(wave)
          wave = wave[srt]
          fx = fx[srt]
          sig = sig[srt]
      endif

      ;; Set vel array
      x_pixminmax, wave, wrest, zabs, vmnx[0], vmnx[1], PIXMIN=pmn, $
        PIXMAX=pmx, VELO=velo

      ;; Y-axis
      if npx EQ 2 AND (ii+1) MOD ntot GT npy then flg_yax = 1 else flg_yax=0

      ;; Plot
      if ii GT npy-1 then ysty=5 else ysty=1
      if (ii NE nlin-1) AND (ii NE npy-1 ) AND $
        ((ii+1) MOD (npx*npy) NE 0) then begin
          spaces = replicate('!17 ',30)
          plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, $
            yrange=ymnx, xtickn=spaces, xmargin=xmrg, $
            ymargin=[0,0], NODATA=nblnd, $
            charsize=csize, psym=10, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=ysty, thick=1.
      endif else begin
          plot, velo[pmn:pmx], fx[pmn:pmx], xrange=vmnx, $
            yrange=ymnx, xmargin=xmrg, ymargin=[0,0], NODATA=nblnd, $
            charsize=csize, psym=10, background=clr.white, color=clr.black, $
            xstyle=1, ystyle=ysty, thick=1., xtickinterval=xtint, $
            xtitle='!17 Relative Velocity (km/s)'
      endelse
      if ii GT npy-1 then begin
          spaces = replicate('!17 ',30)
          axis, yrange=ymnx, ystyle=1, yaxis=1, charsize=csize
          axis, yrange=ymnx, ystyle=1, yaxis=0, ytickn=spaces
      endif

      ;; BLENDS
      nlow = pmn
      for jj=0L,nblnd-1 do begin
          readf, 11, vmin, vmax, FORMAT='(f,f)'
          mn = min(abs(velo-vmin),pixmin)
          mn = min(abs(velo-vmax),pixmax)

          ;; Plot good
          nlow = nlow < pixmin
          oplot, velo[nlow:pixmin], fx[nlow:pixmin], color=clr.black, psym=10, $
            thick=1.
          ;; Plot blend
          oplot, velo[pixmin:pixmax], fx[pixmin:pixmax], color=clr.orange, $
            psym=10, linestyle=2, thick=1.
          ;; End
          if jj EQ nblnd-1 then begin
              if pixmax LT pmx then $
                oplot, velo[pixmax:pmx], fx[pixmax:pmx], color=clr.black, $
                psym=10, thick=1.
          endif
          ;; Set nlow
          nlow = pixmax
      endfor
          

      ;; Labels
      getfnam, wrest, fv, nam
      xyouts, 0.07*(vmnx[1]-vmnx[0])+vmnx[0], $
        ymnx[0]+ (ymnx[1]-ymnx[0])*0.10, $
        strtrim(nam,2), color=clr.black, charsize=LSIZE

      if (ii+1) MOD ntot EQ 1 then begin
          if npx EQ 1 then $
            xyouts, -0.02, 0.5, '!17 Normalized Flux', $
            alignment=0.5, ORIENTATION=90., /normal, charsize=2.0 $
          else $
            xyouts, -0.03, 0.5, '!17 Normalized Flux', $
            alignment=0.5, ORIENTATION=90., /normal, charsize=2.0 
      endif
      
      ;; Lines
      oplot, [0., 0.], ymnx, color=clr.blue, linestyle=2
      oplot, [-10000., 10000.], [0.,0.], color=clr.red, linestyle=3
      oplot, [-10000., 10000.], [1.,1.], color=clr.green, linestyle=3

  endfor

      

  if keyword_set( PSFILE ) then begin
      ps_close, /noprint, /noid
      device, decompose=1
      !p.thick = 1
      !p.charthick = 1
  endif

  !p.multi=[0,1,1]

  return
end
