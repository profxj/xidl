;+ 
; NAME:
; boss_objinf   
;    Version 1.1
;
; PURPOSE:
;  Routine to quickly parse the SDSS QSO database.    
;
; CALLING SEQUENCE:
; boss_objinf, names, SUMMF=, /PLOT, DATDIR=, $
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
;  SUMMF=  -- Summary file [default: 'boss2_qsoi205.fits']
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
;  boss_objinf, '0571-52286-276', /plot
;  boss_objinf, [463, 95], filnm=filnm
;  boss_objinf, ['12:00:21.3', '+32:22:12'], /nmrad
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   16-Nov-2013 Written by JXP
;-
;------------------------------------------------------------------------------

pro boss_objinf, names, SUMMF=SUMMF, PLOT=plot, DATDIR=datdir, $
                 FILNM=sfil, CHKDF=chkdf, PATH=path, ZEM=zem,$
                 RMAG=rmag, RA=ra, DEC=dec, SILENT=silent, NMRAD=nmrad, $
                 IMAG=imag, ZMAG=zmag, NM_QSO=nm_qso, $
                 DR1=dr1, ZIN=zin, PLATE=plate, FIBID=fibid, DR2=dr2, $
                 MJD=mjd, FLGSDSS=flgboss, MOCK=mock, DR3=dr3, LLS=lls, DR4=dr4, $
                 BOSS=boss, OLDZEM=oldzem, DR5=dr5, IDEC=idec, _EXTRA=extra,$

                 DR6=dr6, QSO=qso, PHOTOM=photom, DR7=dr7, OUT_PHOTO=out_photo,$
                 MKCHRT=mkchrt, OLD_SDSS=old_boss, TOLER_RAD=rad_toler

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'boss_objinf, names, /PLOT, CHKDF=, DATDIR=, /NMRAD, FILNm=, RA=, DEC=, /SILENT'
	print, '    ZIN=, SDSS=, FLGSDSS=, '
	print, '    RMAG=, ZEM=, /DR1, PLATE=, FIBID=, /DR2, /DR3, /DR4, /DR5, /DR6, /OLDZEM, /QSO, /PHOTOM, OUT_PHOTO=, TOLER_RAD= [v1.1]'
      return
  endif 

  sfil = ''
  
  if size(names,/type) EQ 7 then NMRAD=1

  if not keyword_set(RAD_TOLER) then RAD_TOLER = 0.005
  if keyword_set(IDEC) then NMRAD = 1
  if not keyword_set( QLIM ) then qlim = 5.

  path = getenv('BOSSPATH')+'/DR10/'
  summf = path+'BOSSLyaDR10_cat_v2.1.fits'
  datdir = path+'BOSSLyaDR10_spectra_v2.1/'
  ;if keyword_set( DR7 ) then begin
  ;    path = getenv('SDSSPATH')+'/DR7_QSO/' 
  ;    summf = path+'dr7_qso.fits'
  ;    datdir = path+'spectro/1d_26/'
  ;endif

  flgboss = 1
  ;;  Summary list
  if not keyword_set( BOSS ) then boss = xmrdfits(summf, 1, /silent)

  ;; Checked file
  if keyword_set( CHKDF ) then schk = xmrdfits(chkdf, 1, /silent)

  ;; Loop
  if size(names,/type) NE 7 AND not keyword_set(NMRAD) then begin
      if n_elements(names) NE 2 then begin
          print, 'boss_objinf: Expecting [plate, fiber] input'
          print, 'boss_objinf: Given -- ', names
          print, 'boss_objfin: Returning..'
          flgboss=0
          return
      endif
      ;; Get mjd
      a = where(boss.plate EQ names[0] AND boss.fiberid EQ names[1],na)
      if na EQ 0 then begin
          if not keyword_set(SILENT) then $
            print, 'boss_objinf: Target '+strtrim(names,2)+ $
                   ' was not found! Returning..'
          flgboss = 0
          return
      endif
      mjd = boss[a[0]].mjd

      ;; Fill up other stuff
      plate = names[0]
      pnm = strtrim(names[0],2)
      if names[0] LT 1000 then pnm='0'+pnm
      fnm = x_padstr(strtrim(names[1],2), 4, '0', /rev)
      ;;
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
              fnm = x_padstr(strmid(names[q],11,3),4,'0',/rev)
              fid = long( fnm )
          endif

          ;; Find the ID
          iobj = where(boss.plate EQ plate AND boss.mjd EQ mjd and $
                       boss.fiberid EQ fid, nobj)
      endif else begin
          ;; JNM
          if flg_nm EQ 2 then begin
              if not keyword_set(IDEC) then $
                x_radec, names[0], names[1], rad, decd $
              else begin
                  rad = names[0]
                  decd = names[1]
              endelse
              iobj = where(abs(boss.RA-rad) LT 0.001 AND $
                           abs(boss.dec-decd) LT 0.001, nobj)
              if nobj eq 0 then begin  
                 ;; Try expanding search
                 iobj = where(abs(boss.RA-rad) LT RAD_TOLER AND $
                              abs(boss.dec-decd) LT RAD_TOLER, nobj)
                 if nobj gt 1 then begin
                    ;; Take closest
                    print,'boss_objinf: Expanded search; find closest ',names[q]
                    dum = min(sqrt((boss[iobj].RA-rad)^2 + $
                                   (boss[iobj].dec-decd)^2),imn)
                    iobj = iobj[imn]
                 endif 
              endif 
          endif
      endelse
      ;;
      case nobj of
          0: begin
              flg = 0
              print, 'boss_objinf: Target ',names[q],' was not found!'
              print, 'boss_objinf: Returning..'
              plate = 0L
              flgboss = 0
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
          plate = boss[iobj].plate
          fibid = boss[iobj].fiberid
          ;; Names
          pnm = strtrim(boss[iobj].plate,2)
          if boss[iobj].plate LT 1000L then pnm = '0'+pnm
          fnm = x_padstr(strtrim(boss[iobj].fiberid,2),4,'0',/rev)
          mjd = boss[iobj].mjd
      endif

      if flg EQ 0 then continue
          

      ;; Print to screen
      x_radec, ra, dec, boss[iobj].RA, boss[iobj].dec, /flip
      if flg_nm EQ 1 then nm = names[q] else nm = pnm+'-'+strtrim(mjd,2)+'-'+fnm
      mag = boss[iobj].psfmag[2] 

      ;; QSO name
      nm_qso = 'J'+strjoin(strsplit(ra,':',/extract))+strjoin(strsplit(dec,':',/extract))

      ;; Redshift
      zem = boss[iobj].z_vi

      if not keyword_set( SILENT ) then $
         print, nm, ra, dec, mag,  zem, $
                FORMAT='(a14,5x,a11,1x,a11,1x,f5.2,1x,f5.3)' 
      
      if arg_present(RMAG) then rmag = mag
      if arg_present(IMAG) then imag = boss[iobj].psfmag[3]
      if arg_present(ZMAG) then zmag = boss[iobj].psfmag[4]

      ;; Chart
      if keyword_set(mkchrt) then begin
         openw,1,'sdss_tmp.inp'
         printf,1,nm+' '+string(sdss[iobj].RA)+' '+string(sdss[iobj].dec)
         close,1
         x_fndchrt,'sdss_tmp.inp',/deci
         spawn,'rm sdss_tmp.inp'
         spawn,'rm tmp.fits'
      endif

      ;; Data file
      sfil = datdir+pnm+'/speclya-'
      sfil = sfil+pnm+'-'+strtrim(mjd,2)+'-'+fnm+'.fits'
      ;; Plot
      if keyword_set( PLOT ) then begin
          if keyword_set( ZIN ) then begin
              if not keyword_set(LLS) then QAL = 1 
          endif
          ;; Plot
          x_specplot, sfil, inflg=19, /block, ZIN=zin, QAL=qal, $
                      LLS=lls, QSO=qso, _EXTRA=extra
      endif
          
  endfor
      
  

return

end
