;+ 
; NAME:
; sdss_objinf   
;    Version 1.1
;
; PURPOSE:
;  Routine to quickly parse the SDSS QSO database.    
;
; CALLING SEQUENCE:
; sdss_objinf, names, SUMMF=, /PLOT, DATDIR=, $
;                FILNM=, CHKDF=, PATH=, ZEM=,$
;                RMAG=, RA=, DEC=, /SILENT, /NMRAD
;
; INPUTS:
;   names -- List of filenames or 2-element array giving plate and fiber
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  SUMMF=  -- Summary file [default: 'sdss2_qsoi205.fits']
;  PATH=   -- Path to SDSS data [default: 'SDSSPATH'+'/DR2_QSO']
;  DATDIR= -- Data directory [default: 'spectro/1d_23']
;  /NMRAD  -- Input is ['RA:::', 'DED:::']
;  /IDEC   -- Input is [rad, decd]
;  /DR1    -- Check the DR1 list (r <= 19.5)
;  ZIN=    -- Input redshift for plotting
;  SDSS=   -- Structure contaning the QSO info
;  /LLS    -- Use the LLS line list
;  /QSO    -- Use the quasar line list
;  /PHOTOM -- Print all photometry
;  RAD_TOLER -- Tolerance for a QSO search
;
; OPTIONAL OUTPUTS:
;  FILNM=  -- File name of the data
;  RA=     -- RA of the QSO
;  DEC=    -- DEC of the QSO
;  FLGSDSS= -- Flag signifying success (1) or failure (0)
;  OUT_PHOTO= -- Array of PSF magnitudes
;  ZEM=    -- Emission redshift
;  NM_QSO  -- QSO name, e.g. J121519.42+424851.0
;
; COMMENTS:
;
; EXAMPLES:
;  sdss_objinf, '0571-52286-276', /plot
;  sdss_objinf, [463, 95], filnm=filnm
;  sdss_objinf, ['12:00:21.3', '+32:22:12'], /nmrad
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   20-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------

pro sdss_objinf, names, SUMMF=SUMMF, PLOT=plot, DATDIR=datdir, $
                 FILNM=sfil, CHKDF=chkdf, PATH=path, ZEM=zem,$
                 RMAG=rmag, RA=ra, DEC=dec, SILENT=silent, NMRAD=nmrad, $
                 IMAG=imag, ZMAG=zmag, NM_QSO=nm_qso, $
                 DR1=dr1, ZIN=zin, PLATE=plate, FIBID=fibid, DR2=dr2, $
                 MJD=mjd, FLGSDSS=flgsdss, MOCK=mock, DR3=dr3, LLS=lls, DR4=dr4, $
                 SDSS=sdss, OLDZEM=oldzem, DR5=dr5, IDEC=idec, _EXTRA=extra,$

                 DR6=dr6, QSO=qso, PHOTOM=photom, DR7=dr7, OUT_PHOTO=out_photo,$
                 MKCHRT=mkchrt, OLD_SDSS=old_sdss, TOLER_RAD=rad_toler

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'sdss_objinf, names, /PLOT, CHKDF=, DATDIR=, /NMRAD, FILNm=, RA=, DEC=, /SILENT'
	print, '    ZIN=, SDSS=, FLGSDSS=, '
	print, '    RMAG=, ZEM=, /DR1, PLATE=, FIBID=, /DR2, /DR3, /DR4, /DR5, /DR6, /OLDZEM, /QSO, /PHOTOM, OUT_PHOTO=, TOLER_RAD= [v1.1]'
      return
  endif 

  sfil = ''
  
  if size(names,/type) EQ 7 then NMRAD=1

  if not keyword_set(RAD_TOLER) then RAD_TOLER = 0.005
  if keyword_set(IDEC) then NMRAD = 1
  if not keyword_set( QLIM ) then qlim = 5.

  ;; Previous DR releases
  if keyword_set( DR1 ) then begin
      path = getenv('SDSSPATH')+'/DR1_QSO/'
      summf = path+'summ_qso_dr1_R195.fit'
      datdir = path+'spectro/1d_20/'
  endif
  if keyword_set( DR2 ) then begin
      path = getenv('SDSSPATH')+'/DR2_QSO/' 
      summf = path+'sdss2_qsoi205.fits'
      datdir = path+'spectro/1d_23/'
  endif
  if keyword_set( DR3 ) then begin
      path = getenv('SDSSPATH')+'/DR3_QSO/' 
      summf = path+'dr3_qso.fits'
      datdir = path+'spectro/1d_23/'
  endif
  if keyword_set( DR4 ) then begin
      path = getenv('SDSSPATH')+'/DR4_QSO/' 
      summf = path+'dr4_qso.fits'
      datdir = path+'spectro/ss_tar_23/'
  endif
  if keyword_set( DR5 ) then begin
      path = getenv('SDSSPATH')+'/DR5_QSO/' 
      summf = path+'dr5_qso.fits'
      datdir = path+'spectro/1d_23/'
  endif
  if keyword_set( DR6 ) then begin
      path = getenv('SDSSPATH')+'/DR6_QSO/' 
      summf = path+'dr6_qso.fits'
      datdir = path+'spectro/1d_25/'
   endif
  if keyword_set( DR7 ) then begin
      path = getenv('SDSSPATH')+'/DR7_QSO/' 
      summf = path+'dr7_qso.fits'
      datdir = path+'spectro/1d_26/'
  endif
  if keyword_set( MOCK ) then begin
      path = getenv('SDSSPATH')+'/Mock/' 
      summf = path+'mock_qso.fits'
      datdir = path+'spectro/1d_23/'
  endif
  if not keyword_set( PATH ) then $
    path = getenv('SDSSPATH')+'/DR7_QSO/' 
  if not keyword_set( SUMMF ) then $
    summf = path+'dr7_qso.fits'
  if not keyword_set( DATDIR ) then $
    datdir = path+'spectro/1d_26/'

  flgsdss = 1
  ;;  Summary list
  if not keyword_set( SDSS ) then sdss = xmrdfits(summf, 1, /silent)

  ;; Checked file
  if keyword_set( CHKDF ) then schk = xmrdfits(chkdf, 1, /silent)

  ;; Loop
  if size(names,/type) NE 7 AND not keyword_set(NMRAD) then begin
      if n_elements(names) NE 2 then begin
          print, 'sdss_objinf: Expecting [plate, fiber] input'
          print, 'sdss_objinf: Given -- ', names
          print, 'sdss_objfin: Returning..'
          flgsdss=0
          return
      endif
      ;; Get mjd
      a = where(sdss.plate EQ names[0] AND sdss.fiberid EQ names[1],na)
      if na EQ 0 then begin
          if not keyword_set(SILENT) then $
            print, 'sdss_objinf: Target '+strtrim(names,2)+ $
                   ' was not found! Returning..'
          flgsdss = 0
          return
      endif
      mjd = sdss[a[0]].mjd

      ;; Fill up other stuff
      plate = names[0]
      pnm = strtrim(names[0],2)
      if names[0] LT 1000 then pnm='0'+pnm
      fnm = strtrim(names[1],2)
      if names[1] LT 10 then fnm = '00'+fnm else $
        if names[1] LT 100 then fnm='0'+fnm
      fid = names[1]
      flg_nm = 0
      nnm = 1
  endif else begin
      if keyword_set( NMRAD ) then begin
          flg_nm = 2
          nnm = 1
      endif else begin
          nnm = n_elements(names)
          flg_nm = 1
      endelse
  endelse

  if not keyword_set( SILENT ) then $
    print, 'Name', 'RA(2000)', 'DEC(2000)', 'Mag(R)', 'zobj', $
    FORMAT='(4x,a4,13x,a8,3x,a9,1x,a6,1x,a4)'
  for q=0L,nnm-1 do begin

      if flg_nm LE 1 then begin
          ;; Parse the name
          if flg_nm EQ 1 then begin
              pnm = strmid(names[q],0,4)
              plate = long( pnm )
              mjd = long( strmid(names[q],5,5) )
              fnm = strmid(names[q],11,3)
              fid = long( fnm )
          endif

          ;; Find the ID
          iobj = where(sdss.plate EQ plate AND sdss.mjd EQ mjd and $
                       sdss.fiberid EQ fid, nobj)
      endif else begin
          ;; JNM
          if flg_nm EQ 2 then begin
              if not keyword_set(IDEC) then $
                x_radec, names[0], names[1], rad, decd $
              else begin
                  rad = names[0]
                  decd = names[1]
              endelse
              iobj = where(abs(sdss.raobj-rad) LT 0.001 AND $
                           abs(sdss.decobj-decd) LT 0.001, nobj)
              if nobj eq 0 then begin  
                 ;; Try expanding search
                 iobj = where(abs(sdss.raobj-rad) LT RAD_TOLER AND $
                              abs(sdss.decobj-decd) LT RAD_TOLER, nobj)
                 if nobj gt 1 then begin
                    ;; Take closest
                    print,'sdss_objinf: Expanded search; find closest ',names[q]
                    dum = min(sqrt((sdss[iobj].raobj-rad)^2 + $
                                   (sdss[iobj].decobj-decd)^2),imn)
                    iobj = iobj[imn]
                 endif 
              endif 
          endif
      endelse
      ;;
      case nobj of
          0: begin
              flg = 0
              print, 'sdss_objinf: Target ',names[q],' was not found!'
              print, 'sdss_objinf: Returning..'
              plate = 0L
              flgsdss = 0
              return
          end
          1: flg=1 
          else: begin
              if keyword_set(NMRAD) then iobj = iobj[0] else stop
              flg = 1
          end
      endcase
      
      if flg_nm EQ 2 then begin
          ;; ID values
          plate = sdss[iobj].plate
          fibid = sdss[iobj].fiberid
          ;; Names
          pnm = strtrim(sdss[iobj].plate,2)
          if sdss[iobj].plate LT 1000L then pnm = '0'+pnm
          fnm = strtrim(sdss[iobj].fiberid,2)
          if sdss[iobj].fiberid LT 10 then fnm = '00'+fnm else $
            if sdss[iobj].fiberid LT 100 then fnm='0'+fnm
          mjd = sdss[iobj].mjd
      endif

      if flg EQ 0 then continue
          

      ;; Print to screen
      x_radec, ra, dec, sdss[iobj].raobj, sdss[iobj].decobj, /flip
      if flg_nm EQ 1 then nm = names[q] else nm = pnm+'-'+strtrim(mjd,2)+'-'+fnm
      if tag_exist(sdss[iobj], 'MAG_R') then mag = sdss[iobj].mag_r $
        else mag = sdss[iobj].psf_r

      ;; QSO name
      nm_qso = 'J'+strjoin(strsplit(ra,':',/extract))+strjoin(strsplit(dec,':',/extract))

      if keyword_set(OLDZEM) then begin
         if not keyword_set(OLD_SDSS) then begin
            path = getenv('SDSSPATH')+'/DR7_QSO/' 
            summf = path+'orig_dr7_qso.fits'
            old_sdss = xmrdfits(summf,1)
         endif
         old_iobj = where(old_sdss.plate EQ plate and $
                      old_sdss.fiberid EQ fid, nni)
         if nni NE 1 then stop
         zem = old_sdss[old_iobj].z 
      endif else zem = sdss_newzem(sdss[iobj].plate, sdss[iobj].fiberid)
      if zem LT 0.001 then zem = sdss[iobj].z

      if not keyword_set( SILENT ) then $
        print, nm, ra, dec, mag,  zem, $
        FORMAT='(a14,5x,a11,1x,a11,1x,f5.2,1x,f5.3)' 

      if arg_present(RMAG) then rmag = mag
      if arg_present(IMAG) then imag = sdss[iobj].psf_i
      if arg_present(ZMAG) then zmag = sdss[iobj].psf_z

      if keyword_set(PHOTOM) or arg_present(OUT_PHOTO) then begin
          print, 'psf_u = ', sdss[iobj].psf_u
          print, 'psf_g = ', sdss[iobj].psf_g
          print, 'psf_r = ', sdss[iobj].psf_r
          print, 'psf_i = ', sdss[iobj].psf_i
          print, 'psf_z = ', sdss[iobj].psf_z
          out_photo = [sdss[iobj].psf_u, $
                       sdss[iobj].psf_g,$
                       sdss[iobj].psf_r, $
                       sdss[iobj].psf_i, $
                       sdss[iobj].psf_z]
       endif

      ;; Checked file
      if keyword_set( CHKDF ) then begin
          a = where(abs(sdss[iobj].raobj - schk.ra) LT 0.0001 AND $
                    abs(sdss[iobj].decobj - schk.dec) LT 0.0001 AND $
                    schk.quality GT qlim, na)
          if na NE 0 then begin
              print, 'zabs', 'f_m', 'f_N', 'NHI', 'qual', $
                FORMAT='(3x,a4,1x,a3,1x,a3,2x,a3,3x,a4)'
              for ii=0L,na-1 do begin
                  print, schk[a[ii]].zabs, schk[a[ii]].flg_mtl, $
                    schk[a[ii]].flg_NHI, $
                    schk[a[ii]].NHI, schk[a[ii]].quality, $
                    FORMAT='(3x,f5.3,1x,i2,1x,i2,2x,f5.2,2x,f4.1)'
              endfor
          endif
      endif
      ;;JMO ADD
      if keyword_set(mkchrt) then begin
         openw,1,'sdss_tmp.inp'
         printf,1,nm+' '+string(sdss[iobj].raobj)+' '+string(sdss[iobj].decobj)
         close,1
         x_fndchrt,'sdss_tmp.inp',/deci
         spawn,'rm sdss_tmp.inp'
         spawn,'rm tmp.fits'
      endif

      ;; Plot
      if keyword_set(DR4) then  sfil = datdir+pnm+'/spSpec/'+'spSpec-' $
      else sfil = datdir+pnm+'/1d/'+'spSpec-' 
      sfil = sfil+strtrim(mjd,2)+'-'+pnm+'-'+fnm+'.fit'
      if keyword_set( PLOT ) then begin
          if keyword_set( ZIN ) then begin
              if not keyword_set(LLS) then QAL = 1 
          endif
          ;; Plot
          x_specplot, sfil, inflg=5, /block, ZIN=zin, QAL=qal, $
                      LLS=lls, QSO=qso, _EXTRA=extra
      endif
          
  endfor
      
  

return

end
