pro mage_arc, arcraw, ordr_str=ordr_str, arcimg=arcimg, $
              INTER=inter, LINLIST=linlist, NOECHO=noecho, $
              CHK=chk, CLOBBER=clobber, SIGREJ=sigrej, DEBUG=debug,$
              PINTER=pinter, IPSIG=ipsig, IFSIG=ifsig, outarc=outarc, $
              GUESSARC=guessarc, SHFTPRM=shftprm, qafil=qafil 


  ;; Optional Keywords
  if not keyword_set( TRCNSIG ) then trcnsig = 3.
  if not keyword_set( RADIUS ) then radius = 2.0
 
  if not keyword_set( LINLIST ) then $
     linlist = getenv('MAGE_DIR')+'/Calib/mm_thar_mod.lst' 
  if not keyword_set( PKWDTH ) then pkwdth = 4L
  if not keyword_set( NOLOG ) then FLOG = 1 else FLOG = 0
  
  ;; Set paramters for peak finding 
  hdr = xheadfits(arcraw)
  slit = FLOAT(SXPAR(hdr, 'SLITNAME'))

  pkwdth=(slit/0.3)*1.2d ;; pixels per slit times 1.2
  TOLER = FLOOR(PKWDTH/3.)
  MXOFF=TOLER
  ;; Lower peak finding sigma threshold for bluest and reddest order first pass
  psig = fltarr(15) + 20.0
  psig[0]=7.0
  psig[1]=7.0
  psig[13]=10.0
  psig[14]=10.0
  ;; final pass
  fsig = fltarr(15) + 10.0
  fsig[0]=5.0
  fsig[1]=5.0
  fsig[13]=7.0
  fsig[14]=5.0
  ;; sigma rejection for wavelength fits
  sigrej=2.0D

  ;; Does the Arc Image already exist?  If so, no need to re-do
  if (keyword_set(arcimg)) then begin
     arcimg_fil = arcimg
  endif else begin
     file_arr = strsplit(arcraw, "mage", /extract)
     arcimg_fil = 'Arcs/ArcImg'+strtrim(file_arr[n_elements(file_arr)-1],2)
     if (x_chkfil(arcimg_fil, /silent) GT 0 AND NOT KEYWORD_SET(CLOBBER)) then begin
        print, "Wavelength file "+arcimg_fil+" already exists. Skipping..."
        return
     endif
     outarc = arcimg_fil
  endelse

  IF NOT FILE_TEST('Arcs', /DIR) THEN FILE_MKDIR, 'Arcs'
  
 ;All the arc files are processed
  mage_proc, arcraw, arcframe, arcinvvar
  mwrfits, arcframe, 'Arcs/Arc'+strtrim(file_arr[n_elements(file_arr)-1],2) $
           , /create
  mwrfits, arcinvvar, 'Arcs/Arc'+strtrim(file_arr[n_elements(file_arr)-1],2)

  if not (keyword_set(GUESSARC)) then $
     guessarc = getenv('MAGE_DIR')+'/Calib/MagE_wvguess_jfh.idl'
  out_fil = 'tmp.idl' ; This is a temporary location for storing the solution

   ;Archived wavelength image...pixel x,y to log lambda
  ;;wvarchive = getenv('MAGE_DIR')+'/Calib/wavelength_new.fits'

  ;;wv = xmrdfits(wvarchive, 0, /fscale)

 
  qafil = 'Arcs/Arc1d_qa'+strtrim(file_arr[n_elements(file_arr)-1],2)+'.ps'
  
                                ;this fits the arc lines down the
                                ;center of each order, and uses an
                                ;archived solution (guessarc) to get
                                ;in the ballpark.  Returns the answer
                                ;to tmp.idl (out_fil)
pinter=0
  result = m_fitarc( 'Arcs/Arc'+strtrim(file_arr[n_elements(file_arr)-1],2) $

                     , ordr_str, out_fil, INTER=0, LINLIST=linlist $ 
                     , CHK=chk, CLOBBER=1, PSIG=psig, FSIG=fsig $
                     , PINTER=pinter, DEBUG=debug, SIGREJ=sigrej $
                     , GUESSARC=guessarc, QAFIL=qafil, NOLOG=nolog $
                     , PKWDTH=PKWDTH, TOLER=TOLER, FIT_MSK=fit_msk $
                     , BCKWD=bckwd, /THIN,/NOEXTRAP)


  ;;Sigma rejection 
  qafil = 'Arcs/Arc2d_qa'+strtrim(file_arr[n_elements(file_arr)-1],2)+'.ps'

  ;;note that nycoeff and nocoeff are
  ;;just going to the defaults here  This
  ;;step fits the arc lines using the
  ;;out_fil from the above step
  sigrej_2d=2.5

  rslt = m_fit2darc( "Arcs/Arc"+strtrim(file_arr[n_elements(file_arr)-1],2) $
                     ,ordr_str, out_fil,nycoeff=6, nocoeff=6 $
                     ,CHKRES=chkres , SIGREJ=sigrej_2d, OUT_STR=out_str $
                     ,CLOBBER=clobber, DEBUG=debug, QAFIL=qafil) ;, /yval)

  mwrfits, out_str, 'Arcs/Arc2d'+strtrim(file_arr[n_elements(file_arr)-1],2)
  qafil = 'Arcs/Trcarc_qa'+strtrim(file_arr[n_elements(file_arr)-1],2)+'.ps'
  out_fil = 'Arcs/TraceArc'+strtrim(file_arr[n_elements(file_arr)-1],2)
  ;;Finds all the lines, traces them and
  ;;that the output is now directed to TraceArc.fits
  rslt = m_echtrcarc('Arcs/Arc'+strtrim(file_arr[n_elements(file_arr)-1],2) $
                     , ordr_str, out_fil, CLOBBER=clobber, INIO=inio $
                     , QAFIL=qafil, szccd=[1024L,2048L])
  arc_fil = 'Arcs/Arc'+strtrim(file_arr[n_elements(file_arr)-1],2)
  trc_fil = out_fil
  qafil = 'Arcs/Trcfit_qa'+strtrim(file_arr[n_elements(file_arr)-1],2)+'.ps'
  out_fil = 'Arcs/Trcfit'+strtrim(file_arr[n_elements(file_arr)-1],2)
  ordr_fil = 'OStr_new_mage.fits'


  ;;Output is now Trcfit.  This routine
  ;;fits the tilt of the lines, and this
  ;;information is used to construct a
  ;;pixel -> wavelength map.
  rslt = m_fittrcarc(arc_fil, trc_fil, ordr_str, out_fil, qafil $
                     ,CHK=chk, CLOBBER=clobber,ORDR_FIL=ordr_fil, _EXTRA=extra) 
  ;;  Final step!
  fil_fittrc = out_fil
  ;; make 2-d wavelength image
  rslt = m_mkaimg(arc_fil, ordr_str $
                  , 'Arcs/Arc2d'+strtrim(file_arr[n_elements(file_arr)-1],2) $
                  , fil_fittrc, arcimg_fil, CHK=chk, CLOBBER=clobber $
                  ,BAD_ANLY=bad_anly) 

end
