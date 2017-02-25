 ;+ 
; NAME:
; mike_subbias   
;     Version 1.2
;
; PURPOSE:
;    Remove bias from all exposures listed in the "mike" structure.
;    The main driver is mike_suboscan which strips the image of the
;    overscan and bias regions (see that program for a full
;    description).  The mike_subbias routine will also remove an
;    archived bias image if requested.
;
; CALLING SEQUENCE:
;   
;  mike_subbias, mike, indx, /usebias, /nobiasrow, /clobber, /ARC,
;  /debug, BADROW=
;
; INPUTS:
;   mike  -  MIKE structure
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
;             were used in mike_mkbias (i.e. /nobiasrow )
;             Presently, this step is not recommended.
;  /ARC -  We found there is significant bleeding into the OV region
;         because of the very bright Arc lines.  Setting /REDARC turns
;         off the fitting in the OV region for the red side.  This may
;         be removed when the neutral density filter is introduced.
;  /CLOBBER = Overwrite existing OV files 
;  /DEBUG -- Turn debug mode on
;  BADROW -- Rows identified as anomolous in the overscan region.
;            Generally the result of an electronics 'hiccup'.
;  /NOSCORR = Do not correct for saturation
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
;   mike_subbias, mike, indx
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
;   mike_subbias_sngl -- Subroutine under mike_subbias.  Uses most of
;                        the keywords described above.  Example:  
;                        rslt = mike_subbias_sngl('Raw/mb0020.fits', side)
;   mike_suboscan
;
; REVISION HISTORY:
;   17-Apr-2003 Written by JXP
;   29-Apr-2003 Modified by GP
;   24-Jun-2003 Modified by JXP (dealt with Blue side)
;
;   18-Jul-2003 Modified by RAB - changed bias subtraction method. 
;               Now as follows:  fit oscan cols (left), subtract
;                                fit oscan rows (top) IF REQUESTED, subtract
;                                subtract bias image  IF REQUESTED
;   19-Aug-2003 Added single file mode
;   26-Jun-2004 Allowed retrieval of image only (no writing)
;   04-Apr-2005 Replaced mike_suboscan with x_suboscan   
;                                
;------------------------------------------------------------------------------

function mike_subbias_sngl, rawfil, side, USEBIAS=usebias, $
                  NOBIASROW = nobiasrow, ARC=arc, CLOBBER=clobber, $
                  DEBUG=debug, IMTYP=imtyp, BADROW=badrow, OVIMG=ovimg, $
                  NOFITS=nofits, SILENT=silent, HEAD=head, NOSCORR=noscorr

  colors = GetColor(/Load, Start=1)

  if  N_params() LT 2  then begin 
      print,'Syntax:  ' + $
        'rslt = mike_subbias_sngl(rawfil, side, /noBIASROW, ' + $
        '/USEBIAS, /ARC, OVROOT=, OVIMG=, /NOFITS, ' + $
        '/CLOBBER, /DEBUG, /SILENT, /NOSCORR [v2.0])'
      return, -1
  endif 
  if keyword_set( NOBIASROW ) then stop
  
; Create OV directory (if necessary)
  a = findfile('OV/..', count=count)
  if count EQ 0 then file_mkdir, 'OV'
  
; Optional Keywords
  if not keyword_set( IMTYP ) then begin
      if keyword_set( ARC ) then imtyp = 'ARC' else imtyp = 'UNK'
  endif

  ;; Check for output
  outfil = mike_getfil('ov_fil', SUBFIL=rawfil, /name, CHKFIL=chkf)
  if CHKF NE 0 and not keyword_set( CLOBBER ) then begin
      print,$
        'mike_subbias_sngl: File ', outfil, $
        ' found. Overscan already subtracted.'
      return, 0
  endif
          
  ;; Open Raw image
  raw = xmrdfits(rawfil, 0, head, /silent, /fscale)
  if not keyword_set( SILENT ) then $
    print, 'mike_subbias: Subtracting the overscan for image ', rawfil
  sz = size(raw, /dimensions)
      
  ;; Set cbin, rbin
  cbin = round(2300. / sz[0] )  ;; Modified for 4x2 binning
  rbin = round(4220. / sz[1] )
      
  ;; JXP :: Am saving the entire science image now
  strcol = 0L 
  fincol = (2048L / cbin) - 1
  strrow = 0L 
  finrow = (4096L / rbin) - 1   ; Last row is usually crummy
          
  if keyword_set (USEBIAS) then begin 
      bias = mike_getfil('bias_fil', SIDE=side, SZ=[2048L/cbin,4096L/rbin], $
                         FIL_NM=bias_fil)
      if not keyword_set( SILENT ) then $
        print, 'mike_subbias: Using BIAS file: ', bias_fil
      bias_size=size(bias,/dimensions)
      if not keyword_set( SILENT ) then $
        print, '   Smoothing bias image.'
      ;; DEBUG
      if keyword_set (debug) then begin
          for ii=10,15 do begin
              plot, findgen(100), bias[ii,100:199]
              smth =  smooth(bias[ii,*],59,/edge_truncate)
              oplot, findgen(100), smth[100:199]
              bias[ii,*] = smth
              oplot, findgen(100), bias[ii,100:199], $
                color=colors.red
          endfor
      endif

      for ii = 0,fincol do begin 
          smth =  smooth(bias[ii,*],59,/edge_truncate)
          bias[ii,*] = smth
      endfor
      if keyword_set (debug) then begin
          xatv, bias, /block
          stop
      endif
              
  endif
          
  ;; Work
  if keyword_set (arc) and side EQ 2 then redarc=1 else redarc=0
      
;  mike_suboscan, raw, head, ovimg, rbin, cbin, imtyp, $ 
;    NoBIASROW= nobiasrow, DEBUG= debug, REDARC=redarc, $
;    SVBAD=badrow, SILENT=silent
  ;; Replaced in V2.0
  x_suboscan, raw, head, ovimg, fincol, RBIN=rbin, CBIN=cbin, $
    IMTYP=imtyp, /BIASROW, DEBUG= debug, SKIPOV=redarc, $
    SVBAD=badrow, SILENT=silent, FINROW=finrow
      
  if keyword_set (USEBIAS) then begin 
      ovimg = temporary(ovimg) - bias
      sxaddpar, head, 'BIAS', 'T', bias_fil
  endif
      
  ovimg = ovimg[strcol:fincol,strrow:finrow]
  imsec = '[' +strtrim(strcol,2)+ ':' +$
    strtrim(fincol,2)+ ',' + $
    strtrim(strrow,2)+ ':'   + $
    strtrim(finrow,2)+ ']'
  sxaddpar, head, 'TRIM', 'T', imsec
    
  ;; Correction for blue side
  if (side EQ 1) AND (sxpar(head,'UT-DATE') GT '2004-05-06') $
                 AND (sxpar(head,'UT-DATE') LT '2005-09-21') $
                 AND NOT keyword_set(NOSCORR) then begin
    correction_factor = exp((((ovimg < 35000.)-14000.) > 0)^2/25000.0^2)
    ovimg = (ovimg * correction_factor) < 50000. 
    sxaddpar, head, 'ADU_CORR', 'T', $
                'A2D factor applied between 14k and 35k DN'
    sxaddpar, head, 'ADU_FORM', 'T', '  = exp((DN-14k)^2/25k^2)'
  endif
   
  ;; Flip 'new' CCD 
  if (side EQ 1) AND (sxpar(head,'UT-DATE') GT '2009-03-15') then begin
      ;; Flipping
      print, 'mike_subbias:  Flipping the blue CCD!'
      ovimg = rotate(ovimg,2)
  endif

  ;; Write out bias-subtracted image
  if not keyword_set(NOFITS) then $
    mwrfits, ovimg, outfil, head, /create, /silent
  
  if not keyword_set( SILENT ) then $
    print, 'mike_subbias_sngl: All done!'
  
  return, 1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_subbias, mike, indx, USEBIAS=usebias, NOBIASROW = nobiasrow, $
                  ARC=arc, CLOBBER=clobber, DEBUG=debug, BADROW=badrow, $
                  HEAD=head, NOSCORR=noscorr
                                ; , VIEW=view

  colors = GetColor(/Load, Start=1)

  if  N_params() LT 2  then begin 
      print,'Syntax:  ' + $
        'mike_subbias, mike, indx, /noBIASROW, /USEBIAS, /ARC, OVROOT=, ' + $
        '/CLOBBER, /DEBUG, BADROW=, /NOSCORR [v1.1]'
;      print,'Recommended:  ' + $
;        'mike_subbias, mike, /noBIASROW, /USEBIAS, /CLOBBER '
      return
  endif 
  
  if keyword_set( NOBIASROW ) then stop
  
  cbin = 0L
  rbin = 0L
  
  nindx = n_elements(indx)
  
  ;; Loop
  for q=0,nindx-1 do begin
      side = mike[indx[q]].side
      ;; Check for output
      outfil = mike_getfil('ov_fil', SUBFIL=mike[indx[q]].img_root, $
                           /name, CHKFIL=chkf)
      if CHKF NE 0 and not keyword_set( CLOBBER ) then begin
          print,$
            'mike_subbias: File ', outfil, $
            ' found. Overscan already subtracted.'
          mike[indx[q]].img_ov = outfil
          mike[indx[q]].flg_ov = 1
          continue
      endif
          
      ;; Open Raw image
      raw = xmrdfits(mike[indx[q]].rootpth+mike[indx[q]].img_root, 0, head, $
                     /silent, /fscale)
      print, 'mike_subbias: Subtracting the overscan for image ', $
        mike[indx[q]].rootpth+mike[indx[q]].img_root,  ' ', mike[indx[q]].type
      
      if (cbin ne mike[indx[q]].colbin or rbin ne mike[indx[q]].rowbin) then begin
          cbin = mike[indx[q]].colbin 
          rbin = mike[indx[q]].rowbin
          
          ;; JXP :: Am saving the entire science image now
          strcol = 0L 
          fincol = (2048L / cbin) - 1
          strrow = 0L 
          finrow = (4096L / rbin) - 1 ; Last row is usually crummy
          
          if keyword_set (USEBIAS) then begin 
              bias = mike_getfil('bias_fil', SIDE=mike[indx[q]].side, $
                                 SZ=[2048L/cbin,4096L/rbin], FIL_NM=bias_fil)
              print, 'mike_subbias: Using BIAS file: ', bias_fil
              bias_size=size(bias,/dimensions)
              print, '   Smoothing bias image.'
              ;; DEBUG
              if keyword_set (debug) then begin
                  for ii=10,15 do begin
                      plot, findgen(100), bias[ii,100:199]
                      smth =  smooth(bias[ii,*],59,/edge_truncate)
                      oplot, findgen(100), smth[100:199]
                      bias[ii,*] = smth
                      oplot, findgen(100), bias[ii,100:199], $
                        color=colors.red
                  endfor
              endif


              for ii = 0,fincol do begin 
                  smth =  smooth(bias[ii,*],59,/edge_truncate)
                  bias[ii,*] = smth
              endfor
              
              if keyword_set (debug) then begin
                  xatv, bias, /block
                  stop
              endif
              
          endif
          
      endif

      ;; Work
      if keyword_set (arc) AND mike[indx[q]].type eq 'ARC' $
        and mike[indx[q]].side EQ 2 then redarc=1 else redarc=0
      
;      mike_suboscan, raw, head, ovimg, $
;        mike[indx[q]].rowbin, mike[indx[q]].colbin, mike[indx[q]].type, $ 
;        NoBIASROW= nobiasrow, DEBUG= debug, REDARC=redarc, SVBAD=badrow
      ;; Replaced in V2.0
      x_suboscan, raw, head, ovimg, fincol, RBIN=mike[indx[q]].rowbin, $
        CBIN=mike[indx[q]].colbin, IMTYP=mike[indx[q]].type, $
        /BIASROW, DEBUG= debug, SKIPOV=redarc, $
        SVBAD=badrow, SILENT=silent, FINROW=finrow
      
      if keyword_set (USEBIAS) then begin 
          ovimg = temporary(ovimg) - bias
          mike[indx[q]].flg_anly = 3
          sxaddpar, head, 'BIAS', 'T', bias_fil
      endif
      
      ovimg = ovimg[strcol:fincol,strrow:finrow]
      imsec = '[' +strtrim(strcol,2)+ ':' +$
        strtrim(fincol,2)+ ',' + $
        strtrim(strrow,2)+ ':'   + $
        strtrim(finrow,2)+ ']'
      sxaddpar, head, 'TRIM', 'T', imsec

      ;; Correct the counts
      if (side EQ 1) AND (sxpar(head,'UT-DATE') GT '2004-05-06') $
                     AND (sxpar(head,'UT-DATE') LT '2005-09-21') then begin
        correction_factor = exp((((ovimg < 35000.)-14000.) > 0)^2/25000.0^2)
        ovimg = (ovimg * correction_factor) < 50000. 
        sxaddpar, head, 'ADU_CORR', 'T', $
                'A2D factor applied between 14k and 35k DN'
        sxaddpar, head, 'ADU_FORM', 'T', '  = exp((DN-14k)^2/25k^2) 
      endif
   
      
      mike[indx[q]].img_ov = outfil
      mike[indx[q]].flg_ov = 1

      ;; Flip if the 'new' CCD
      if (side EQ 1) AND (sxpar(head,'UT-DATE') GT '2009-03-15') then begin
          print, 'mike_subbias:  Flipping the blue CCD!'
          ovimg = rotate(ovimg,2)
      endif

      ;; Write out bias-subtracted image
      mwrfits, ovimg, outfil, head, /create, /silent
      
  endfor

  print, 'mike_subbias: All done!'

  return
end
