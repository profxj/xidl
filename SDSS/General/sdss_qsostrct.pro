;+ 
; NAME:
; sdss_qsostrct
;    Version 1.1
;
; PURPOSE:
;  Routine to quickly grab the SDSS QSO structure
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
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  qso = sdss_qsostrct()
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   15-Mar-2010 Written by JXP
;-
;------------------------------------------------------------------------------

function sdss_qsostrct, DR

  
  if not keyword_set(DR) then DR = 'DR7'

  ;; Previous DR releases
  case DR of
      'DR1': begin
          path = getenv('SDSSPATH')+'/DR1_QSO/'
          summf = path+'summ_qso_dr1_R195.fit'
          datdir = path+'spectro/1d_20/'
      end
      'DR2': begin
          path = getenv('SDSSPATH')+'/DR2_QSO/' 
          summf = path+'sdss2_qsoi205.fits'
          datdir = path+'spectro/1d_23/'
      end
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
          endif
      endelse
      ;;
      case nobj of
          0: begin
              flg = 0
              print, 'sdss_objinf: Target '+names[q]+' was not found!'
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

      if keyword_set(OLDZEM) then zem = sdss[iobj].z $
      else zem = sdss_newzem(sdss[iobj].plate, sdss[iobj].fiberid)
      if zem LT 0.001 then zem = sdss[iobj].z

      if not keyword_set( SILENT ) then $
        print, nm, ra, dec, mag,  zem, $
        FORMAT='(a14,5x,a11,1x,a11,1x,f5.2,1x,f5.3)' 

      if arg_present(RMAG) then rmag = mag
      if arg_present(IMAG) then imag = sdss[iobj].psf_i
      if arg_present(ZMAG) then zmag = sdss[iobj].psf_z

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
