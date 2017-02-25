;+ 
; NAME:
; uves_subbias   
;     Version 1.2
;
; PURPOSE:
;    Remove bias from all exposures listed in the "uves" structure.
;    The main driver is uves_suboscan which strips the image of the
;    overscan region (see that program for a full
;    description).  The uves_subbias routine will also remove an
;    archived bias image if requested.
;
; CALLING SEQUENCE:
;   
;  uves_subbias, uves, indx, /usebias, /nobiasrow, /clobber, /ARC,
;  /debug, BADROW=
;
; INPUTS:
;   uves  -  MIKE structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /NOBIASROW= if set, bias row is not used (normally, you should not use it)
;  /USEBIAS = use the bias image in addition to subtracting a fit to the
;             overscan columns (at right).  The code first generates a
;             smoothed version of the bias image.
;             If you do this, you want to be sure the same options
;             were used in uves_mkbias (i.e. /nobiasrow )
;             Presently, this step is not recommended.
;  /CLOBBER = Overwrite existing OV files 
;  /DEBUG -- Turn debug mode on
;  BADROW -- Rows identified as anomolous in the overscan region.
;            Generally the result of an electronics 'hiccup'.
;
; 
; OUTPUTS TO SNGL:
;   /NOFITS -- Do not write a fits file
;   OVIMG=    -- Bias subtracted image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
;
; EXAMPLES:
;   uves_subbias, uves, indx, /usebias, /nobiasrow
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
;   uves_subbias_sngl -- Subroutine under uves_subbias.  Uses most of
;                        the keywords described above.  Example:  
;                        rslt = uves_subbias_sngl('Raw/mb0020.fits',
;                        chip)
;   uves_suboscan
;
; REVISION HISTORY:
;   01-Feb-2005 Written by JXP  (taken from mike_subbias)
;                                
;------------------------------------------------------------------------------

function uves_subbias_sngl, rawfil, side, $
  USEBIAS=usebias, NOBIASROW = nobiasrow, EXTEN=exten, $
  CLOBBER=clobber, DEBUG=debug, IMTYP=imtyp, $
  BADROW=badrow, OVIMG=ovimg, NOFITS=nofits, $
  SILENT=silent, CBIN=cbin, RBIN=rbin, POSCAN=poscan
  

  if  N_params() LT 2  then begin 
      print,'Syntax:  ' + $
        'rslt = uves_subbias_sngl(rawfil, side, /noBIASROW, ' + $
        '/USEBIAS, OVROOT=, OVIMG=, /NOFITS, ' + $
        '/CLOBBER, /DEBUG, /SILENT [v1.2])'
      return, -1
  endif 

  if not keyword_set(EXTEN) then exten = 0L
  ;; Root name (assumes starting with Raw/)
  pos = strpos(rawfil, '.gz')
  if pos LE 0 then pos = strlen(rawfil)
  img_root = strmid(rawfil,4,pos-4)

  outfil = uves_getfil('ov_fil', OBJN=img_root, /name, CHKFIL=chkf)
  if CHKF NE 0 and not keyword_set( CLOBBER ) then begin
      print,$
        'uves_subbias: File ', outfil, $
        ' found. Overscan already subtracted.'
      return, 0
  endif
          
  ;; Open Raw image
  raw = xmrdfits(rawfil, EXTEN, head,  /silent, /fscale)
  print, 'uves_subbias: Subtracting the overscan for image ', rawfil
      
  NAXIS1 = sxpar(head, 'NAXIS1')
  NAXIS2 = sxpar(head, 'NAXIS2')

  if EXTEN GT 0 then begin
      wcen = uves_getwcen(rawfil)
      sxaddpar, head, 'WLEN2', wcen
  endif
  
;  wcen = uves_getwcen(raw_fil)
;  if wcen LT 400 then begin
;      typ = uves_headid(head, 1, SARA=img_root)
;      colbin= sxpar(head,'BINX')
;      rowbin= sxpar(head,'BINY')
;      poscan= [sxpar(head,'PRESCAN'), sxpar(head,'OSCAN')] 
;  endif

  ;; FINCOL
  dats = [ poscan[0], NAXIS1-poscan[1]-1, 0L, NAXIS2-1]
  fincol = dats[1]

  ;; TRIM Vignetting
  if side EQ 1 then $
    dats[3] = dats[3] < (round(2940./rbin) - 1.) $
  else dats[3] = dats[3] < (round(4060./rbin) - 1.) 

  ;; Subtract overscan
  x_suboscan, raw, head, ovimg, fincol+5, $
              DEBUG= debug, SVBAD=badrow, CBIN=cbin,  RBIN=rbin
  
  ;; Trim
;      dats = dats-1
  ovimg = ovimg[dats[0]:dats[1],dats[2]:dats[3]]
;  sxaddpar, head0, 'TRIM', 'T'
;      sxaddpar, head0, 'ODSEC', sdats
  
  ;; Header
  mkhdr, main_head, ovimg
  sxdelpar, main_head, 'END'
  sxdelpar, head, 'NAXIS'
  sxdelpar, head, 'NAXIS1'
  sxdelpar, head, 'NAXIS2'
  sxdelpar, head, 'BITPIX'
  sxdelpar, head, 'BZERO'
  sxdelpar, head, 'SIMPLE'
  sxdelpar, head, 'DATSUM'
  sxdelpar, head, 'DATE'
  
  nhd = n_elements(head)
  for qq=0L,nhd-1 do begin
      if strlen(strtrim(head[qq],2)) GT 0 then $
        main_head = [main_head, head[qq]]
  endfor
  
  mwrfits, ovimg, outfil, main_head, /create, /silent
  
  if not keyword_set( SILENT ) then $
    print, 'uves_subbias_sngl: All done!'

  return, 1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro uves_subbias, uves, indx, USEBIAS=usebias, NOBIASROW = nobiasrow, $
                  CLOBBER=clobber, DEBUG=debug, BADROW=badrow
              ; , VIEW=view

  colors = GetColor(/Load, Start=1)

  if  N_params() LT 2  then begin 
      print,'Syntax:  ' + $
        'uves_subbias, uves, indx, /noBIASROW, /USEBIAS, OVROOT=, ' + $
        '/CLOBBER, /DEBUG, BADROW= [v1.1]'
;      print,'Recommended:  ' + $
;        'uves_subbias, uves, /noBIASROW, /USEBIAS, /CLOBBER '
      return
  endif 
  
  
  cbin = 0L
  rbin = 0L
  
  nindx = n_elements(indx)
  
  ;; Loop
  for q=0,nindx-1 do begin
;      chip = uves[indx[q]].chip
      ;; Check for output
      outfil = uves_getfil('ov_fil', $
                   OBJN=uves[indx[q]].img_root, /name, CHKFIL=chkf)
      if CHKF NE 0 and not keyword_set( CLOBBER ) then begin
          print,$
            'uves_subbias: File ', outfil, $
            ' found. Overscan already subtracted.'
          uves[indx[q]].img_ov = outfil
          uves[indx[q]].flg_ov = 1
          continue
      endif
          
      ;; Open Raw image
      raw = xmrdfits(uves[indx[q]].rootpth+uves[indx[q]].img_root, $
                     uves[indx[q]].exten, head, $
                     /silent, /fscale)
      head0 = xheadfits(uves[indx[q]].rootpth+uves[indx[q]].img_root)
      print, 'uves_subbias: Subtracting the overscan for image ', $
        uves[indx[q]].rootpth+uves[indx[q]].img_root,  ' ', uves[indx[q]].type

      NAXIS1 = sxpar(head, 'NAXIS1')
      NAXIS2 = sxpar(head, 'NAXIS2')

      ;; FINCOL
      dats = [ uves[indx[q]].poscan[0], $
               NAXIS1-uves[indx[q]].poscan[1]-1, 0L, NAXIS2-1]
;      sdats = sxpar(head, 'DATASEC')
      fincol = dats[1]

      ;; TRIM Vignetting
      if uves[indx[q]].side EQ 1 then $
        dats[3] = dats[3] < (round(2940./uves[indx[q]].rowbin) - 1.) $
      else dats[3] = dats[3] < (round(4060./uves[indx[q]].rowbin) - 1.) 

      ;; Subtract overscan
      x_suboscan, raw, head, ovimg, fincol+5, IMTYPE=uves[indx[q]].type, $ 
        DEBUG= debug, SVBAD=badrow, CBIN=uves[indx[q]].colbin, $
        RBIN=uves[indx[q]].rowbin
      
      ;; Trim
;      dats = dats-1
      ovimg = ovimg[dats[0]:dats[1],dats[2]:dats[3]]
      sxaddpar, head0, 'TRIM', 'T'
;      sxaddpar, head0, 'ODSEC', sdats
      
      ;; Write out bias-subtracted image
      uves[indx[q]].img_ov = outfil
      uves[indx[q]].flg_ov = 1

      ;; Header
      mkhdr, main_head, ovimg
      sxdelpar, main_head, 'END'
      sxdelpar, head, 'NAXIS'
      sxdelpar, head, 'NAXIS1'
      sxdelpar, head, 'NAXIS2'
      sxdelpar, head, 'BITPIX'
      sxdelpar, head, 'BZERO'
      sxdelpar, head, 'SIMPLE'
      sxdelpar, head, 'DATSUM'
      sxdelpar, head, 'DATE'

      nhd = n_elements(head)
      for qq=0L,nhd-1 do begin
          if strlen(strtrim(head[qq],2)) GT 0 then $
            main_head = [main_head, head[qq]]
      endfor
      
      mwrfits, ovimg, outfil, main_head, /create, /silent
      
  endfor

  print, 'uves_subbias: All done!'

  return
end
