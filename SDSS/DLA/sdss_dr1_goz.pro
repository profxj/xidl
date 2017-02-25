;+ 
; NAME:
; sdss_dr1_goz
;
; PURPOSE:
;  g(z) prescription for DR1.  OBSOLETE
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
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   27-Feb-2004 Written by SHF
;-
;------------------------------------------------------------------------------
pro sdss_dr1_goz, ntrials=ntrials, QALFIL=qalfil, OUTFIL=outfil, XZFIL=xzfil

  if not keyword_set( XZFIL ) then $
    xzfil = '/u/xavier/SDSS/DR1_QSO/xz_val.fits'
  if not keyword_set( QALFIL ) then $
    qalfil = '/u/xavier/SDSS/DR1_QSO/sdss_DR1_QAL.fits'
  if not keyword_set( OFFL ) then offl = 0.  ;; Offset to start wave

 ;;read in spectra
  x=xmrdfits(qalfil,1,/silent)
  nspectra=n_elements(x)
  f=findgen(5000)*.001
  g=lonarr(5000)


  ;; XZFIL
  xzz = xmrdfits(xzfil, 0, /silent)
  xzv = xmrdfits(xzfil, 0, /silent)

  ;; RA
  ra=fltarr(nspectra)
  dec=fltarr(nspectra)
  zem=fltarr(nspectra)
  if not keyword_set(ntrials) then ntrials=nspectra

  ;;loop over spectra 
  gdqso = [-1L]
  allz1 = [0.d]
  allz2 = [0.d]
  dX = [0.d]
  svra = [0.d]
  svbal = [0]
  svdec = [0.d]
  svzem = [0.d]
  svmag = [0.d]
  svs2n = [0.]
  for i=0,ntrials-1 do begin
      if x[i].z_qso LE 2. then continue

      name=x[i].file_name
;      d1=xmrdfits(name, 0, head, /silent)
      
      ra[i]=x[i].ra
      dec[i]=x[i].dec
      zem[i]=x[i].z_qso
      
      ;; Check for double
      if i gt 0 then begin   
          dbl=where(abs(ra[i]-ra[0:i-1]) le 1E-3 $
                    and abs(dec[i]-dec[0:i-1]) le 1E-3 and $
                    abs(zem[i]-zem[0:i-1]) le .1,ndbl) 
          if ndbl NE 0 then begin
              print, 'sdss_dr1_goz: Repeat quasar!'
          endif
      endif else ndbl=0

      if ndbl le 0 and x[i].start_wave GT 0. then begin
          ;; BAL
          flg = 0
          zstrt = (x[i].start_wave)/1215.6701 -1
          case x[i].flg_bal of
              0: begin
                  z1 = zstrt*1.005 + 0.005  ; Add 1500 km/s to start wave
                  z2 = .99*x[i].z_qso - 0.01  ; 3000 km/s blue of Lya
                  gdi= where(f ge z1 and f le z2, ngd)
                  if ngd ne 0 then g[gdi]=g[gdi]+1
              end
              1: begin
                  z1=(((1060.*(1+x[i].z_qso))/1215.6701)-1) $
                    > (zstrt*1.005 + 0.005)
                  z2=x[i].z_qso - 0.08226  ; Equivalent to 100 Ang
                  gdi=where(f ge z1 and f le z2, ngd)
                  
                  if ngd ne 0 then g[gdi]=g[gdi]+1
              end
              2: begin  ;; Rejected BAL
                  z1 = -1.
                  z2 = -1.
                  flg = -1
              end
              else: flg = -1
          endcase
          if flg NE -1 then begin
              ;; Min on z1 of 2.2 for SDSS
              z1 = z1 > 2.2
              ;;
              gdqso = [gdqso,i]
              allz1 = [allz1, z1]
              allz2 = [allz2, z2]
              svbal = [svbal, x[i].flg_bal]
              svra = [svra, ra[i]]
              svdec = [svdec, dec[i]]
              svzem = [svzem,zem[i]]

              ;; Magnitude and s2n
              svmag = [svmag,x[i].qso_mag]
              svs2n = [svs2n, x[i].snr]

              ;; Calculate dX
              mn = min(abs(xzz-z1),imn1)
              mn = min(abs(xzz-z2),imn2)
              dX = [dx, xzv[imn2]-xzv[imn1]]
          endif
      endif
  endfor

  ;;plot
  x_splot,f,g,/block
  maxg=max(g)
;stop
  old_device=!d.name
  set_plot,'ps'
  device,filename='DR1_goz.ps'
;device,xsize=7,ysize=7,/inches
  plot,f,g,yrange=[0,maxg+1]
  device,/close_file
  set_plot,old_device

  ;; OUT
  if keyword_set( OUTFIL ) then begin
      strct = { $
                gzz: f, $
                gzv: g, $
                ra: svra[1:*], $
                dec: svdec[1:*], $
                flg_bal: svbal[1:*], $
                iqso: gdqso[1:*], $
                mag: svmag[1:*], $
                s2n: svs2n[1:*], $
                z1: allz1[1:*], $
                z2: allz2[1:*], $
                zem: svzem[1:*], $
                dX: dX[1:*] $
              }
      mwrfits,strct,outfil,/create
  endif
  print, 'sdss_dr1_goz: All done!'

end

