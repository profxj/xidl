;+ 
; NAME:
; hires_subbias   
;     Version 1.2
;
; PURPOSE:
;    Remove bias from all exposures listed in the "hires" structure.
;    The main driver is hires_suboscan which strips the image of the
;    overscan region (see that program for a full
;    description).  The hires_subbias routine will also remove an
;    archived bias image if requested.
;
; CALLING SEQUENCE:
;   
;  hires_subbias, hires, indx, /usebias, /nobiasrow, /clobber, /ARC,
;  /debug, BADROW=, RBIN=, CBIN=
;
; INPUTS:
;   hires  -  MIKE structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  CBIN=  -- Column binning
;  RBIN=  -- Row binning
;  /NOBIASROW= if set, bias row is not used (normally, you should not use it)
;  /USEBIAS = use the bias image in addition to subtracting a fit to the
;             overscan columns (at right).  The code first generates a
;             smoothed version of the bias image.
;             If you do this, you want to be sure the same options
;             were used in hires_mkbias (i.e. /nobiasrow )
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
;   hires_subbias, hires, indx, /usebias, /nobiasrow
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
;   hires_subbias_sngl -- Subroutine under hires_subbias.  Uses most of
;                        the keywords described above.  Example:  
;                        rslt = hires_subbias_sngl('Raw/mb0020.fits',
;                        chip)
;   hires_suboscan
;
; REVISION HISTORY:
;   01-Feb-2005 Written by JXP  (taken from mike_subbias)
;                                
;------------------------------------------------------------------------------

function hires_subbias_sngl, rawfil, chip, exten, $
                             USEBIAS=usebias, NOBIASROW = nobiasrow, $
                             CLOBBER=clobber, DEBUG=debug, IMTYP=imtyp, $
                             BADROW=badrow, OVIMG=ovimg, NOFITS=nofits, $
                             SILENT=silent, CBIN=cbin, RBIN=rbin, FRAME=frame

  colors = GetColor(/Load, Start=1)

  if  N_params() LT 3  then begin 
      print,'Syntax:  ' + $
        'rslt = hires_subbias_sngl(rawfil, chip, exten, /noBIASROW, ' + $
        '/USEBIAS, OVROOT=, OVIMG=, /NOFITS, ' + $
        '/CLOBBER, /DEBUG, /SILENT [v1.2])'
      return, -1
  endif 
  
; Create OV directory (if necessary)
  a = findfile('OV/..', count=count)
  if count EQ 0 then file_mkdir, 'OV'
  
; Optional Keywords
  if not keyword_set( IMTYP ) then begin
      if keyword_set( ARC ) then imtyp = 'ARC' else imtyp = 'UNK'
  endif

  if not keyword_set( FRAME ) then begin
      pos1 = strpos(rawfil, 'HI.')
      if pos1 NE -1 then stop
      pos = strpos(rawfil, '.fits')
      frame = long(strmid(rawfil,pos-4,4))
  endif

  ;; Check for output
  outfil = hires_getfil('ov_fil', CHIP=chip, FRAME=frame, /name, CHKFIL=chkf)
  if CHKF NE 0 and not keyword_set( CLOBBER ) then begin
      print,$
        'hires_subbias_sngl: File ', outfil, $
        ' found. Overscan already subtracted.'
      return, 0
  endif
          
  ;; Open Raw image
  raw = xmrdfits(rawfil, exten, head, /silent, /fscale)
  head0 = xheadfits(rawfil)
  if not keyword_set( SILENT ) then $
    print, 'hires_subbias: Subtracting the overscan for image ', rawfil
  sz = size(raw, /dimensions)
      
  ;; Set cbin, rbin
  if not keyword_set( CBIN ) then cbin = round(2048. / sz[0] )
  if not keyword_set( RBIN ) then begin 
     if CHIP NE -1 then rbin = round(4096. / sz[1] ) else $
     rbin = round(2048. / sz[1] ) 
  endif
      
  ;; JXP :: Am saving the entire science image now
  if strmatch(strmid(sxpar(head0,'DETECTOR'),0,20),'HIRES Science Mosaic') then begin
     dats = xregtovec(sxpar(head, 'DATASEC'))
     sdats = sxpar(head, 'DATASEC')
  endif else begin
     pre = sxpar(head, 'PREPIX')
     post = sxpar(head, 'POSTPIX')
     tmp = strsplit(sxpar(head, 'WINDOW'),',',/extract)
     dats = lonarr(4)
     dats[0] = pre
     dats[1] = pre + tmp[3] - 1
     dats[2] = 0
     dats[3] = tmp[4]-1
     sdats = '['+strtrim(dats[0],2)+':'+$
             strtrim(dats[1],2)+','+$
             strtrim(dats[2],2)+':'+$
             strtrim(dats[3],2)+']'
  endelse
  fincol = dats[1]
 
  ;; TRIM Vignetting
  dats[3] = dats[3] < round(3990./rbin)


  ;; TRIM Red chip
  if chip EQ 3L then dats[0] = dats[0] > round(110./cbin)
          
  if keyword_set (USEBIAS) then begin 
      bias = hires_getfil('bias_fil', CHIP=chip, SZ=[2048L/cbin,4096L/rbin], $
                         FIL_NM=bias_fil)
      if not keyword_set( SILENT ) then $
        print, 'hires_subbias: Using BIAS file: ', bias_fil
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

      for ii = 0,fincol-1 do begin 
          smth =  smooth(bias[ii,*],59,/edge_truncate)
          bias[ii,*] = smth
      endfor
      if keyword_set (debug) then begin
          xatv, bias, /block
          stop
      endif
              
  endif
          
  ;; Work
      
  if keyword_set (USEBIAS) then begin 
      ovimg = temporary(ovimg) - bias
      sxaddpar, head, 'BIAS', 'T', bias_fil
  endif

  ;; Ovesrcan
  x_suboscan, raw, head, ovimg, fincol+5, IMTYPE=imtyp, $
    DEBUG= debug, SVBAD=badrow, CBIN=cbin, RBIN=rbin
      
  ;; Trim
  if chip NE -1 then dats = dats-1
  ovimg = ovimg[dats[0]:dats[1],dats[2]:dats[3]]
  sxaddpar, head0, 'TRIM', 'T'
  sxaddpar, head0, 'ODSEC', sdats
      
  ;; Transpose original chip
  if chip EQ -1 then ovimg = rotate(ovimg, 1)

  ;; Header
  mkhdr, main_head, ovimg
  sxdelpar, main_head, 'END'
  sxdelpar, head0, 'NAXIS'
  sxdelpar, head0, 'BITPIX'
  sxdelpar, head0, 'BZERO'
  sxdelpar, head0, 'SIMPLE'
  sxdelpar, head0, 'DATSUM'
  sxdelpar, head0, 'DATE'

  ;; Single chip
  if CHIP EQ -1 then begin
     sxdelpar, head0, 'NAXIS'
     sxdelpar, head0, 'NAXIS1'
     sxdelpar, head0, 'NAXIS2'
     sxdelpar, head0, 'BITPIX'
  endif

  nhd = n_elements(head0)
  for qq=0L,nhd-1 do begin
      if strlen(strtrim(head0[qq],2)) GT 0 then $
        main_head = [main_head, head0[qq]]
  endfor
      
  ;; Write out bias-subtracted image
  if not keyword_set(NOFITS) then $
    mwrfits, ovimg, outfil, main_head, /create, /silent
      
  if not keyword_set( SILENT ) then $
    print, 'hires_subbias_sngl: All done!'

  return, 1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_subbias, hires, indx, USEBIAS=usebias, NOBIASROW = nobiasrow, $
                  CLOBBER=clobber, DEBUG=debug, BADROW=badrow
              ; , VIEW=view

  colors = GetColor(/Load, Start=1)

  if  N_params() LT 2  then begin 
      print,'Syntax:  ' + $
        'hires_subbias, hires, indx, /noBIASROW, /USEBIAS, OVROOT=, ' + $
        '/CLOBBER, /DEBUG, BADROW= [v1.1]'
;      print,'Recommended:  ' + $
;        'hires_subbias, hires, /noBIASROW, /USEBIAS, /CLOBBER '
      return
  endif 
  
  
  cbin = 0L
  rbin = 0L
  
  nindx = n_elements(indx)
  
  ;; Loop
  for q=0,nindx-1 do begin
      chip = hires[indx[q]].chip
      ;; Check for output
      outfil = hires_getfil('ov_fil', $
                   CHIP=hires[indx[q]].chip, $
                   FRAME=hires[indx[q]].frame, /name, CHKFIL=chkf)
      if CHKF NE 0 and not keyword_set( CLOBBER ) then begin
          print,$
            'hires_subbias: File ', outfil, $
            ' found. Overscan already subtracted.'
          hires[indx[q]].img_ov = outfil
          hires[indx[q]].flg_ov = 1
          continue
      endif
          
      ;; Open Raw image
      raw = xmrdfits(hires[indx[q]].rootpth+hires[indx[q]].img_root, $
                     hires[indx[q]].exten, head, $
                     /silent, /fscale)
      head0 = xheadfits(hires[indx[q]].rootpth+hires[indx[q]].img_root)
      print, 'hires_subbias: Subtracting the overscan for image ', $
        hires[indx[q]].rootpth+hires[indx[q]].img_root,  ' ', hires[indx[q]].type

      ;; FINCOL
      if ~ strmatch(hires[0].ccd, 'SINGLE') then begin
         dats = xregtovec(sxpar(head, 'DATASEC'))
         sdats = sxpar(head, 'DATASEC')
      endif else begin
         pre = sxpar(head, 'PREPIX')
         post = sxpar(head, 'POSTPIX')
         tmp = strsplit(sxpar(head, 'WINDOW'),',',/extract)
         dats = lonarr(4)
         dats[0] = pre
         dats[1] = pre + tmp[3] - 1
         dats[2] = 0
         dats[3] = tmp[4]-1
         sdats = '['+strtrim(dats[0],2)+':'+$
                 strtrim(dats[1],2)+','+$
                 strtrim(dats[2],2)+':'+$
                 strtrim(dats[3],2)+']'
      endelse
      fincol = dats[1]

      ;; TRIM Vignetting
      dats[3] = dats[3] < round(3990./hires[indx[q]].rowbin)

      ;; TRIM Red chip
      if chip EQ 3L then dats[0] = dats[0] > round(110./hires[indx[q]].colbin)

      if keyword_set (USEBIAS) then begin 
          bias = hires_getfil('bias_fil', CHIP=hires[indx[q]].chip, $
                              SZ=[2048L/cbin,4096L/rbin], FIL_NM=bias_fil)
          print, 'hires_subbias: Using BIAS file: ', bias_fil
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
          
          
          if keyword_set(SMOOTH) then begin
              for ii = 0,fincol do begin 
                  smth =  smooth(bias[ii,*],59,/edge_truncate)
                  bias[ii,*] = smth
              endfor
          endif
          
          if keyword_set (debug) then begin
              xatv, bias, /block
              stop
          endif
          
      endif
      

      
      if keyword_set (USEBIAS) then begin 
          ovimg = temporary(ovimg) - bias
          hires[indx[q]].flg_anly = 3
          sxaddpar, head0, 'BIAS', 'T', bias_fil
      endif

      ;; Subtract overscan
      x_suboscan, raw, head, ovimg, fincol+5, IMTYPE=hires[indx[q]].type, $ 
        DEBUG= debug, SVBAD=badrow, CBIN=hires[indx[q]].colbin, $
        RBIN=hires[indx[q]].rowbin
      
      ;; Trim
      if ~strmatch(hires[0].ccd,'SINGLE') then dats = dats-1
      ovimg = ovimg[dats[0]:dats[1],dats[2]:dats[3]]
      sxaddpar, head0, 'TRIM', 'T'
      sxaddpar, head0, 'ODSEC', sdats

      ;; Transpose original chip
      if strmatch(hires[0].ccd,'SINGLE') then ovimg = rotate(ovimg, 1)
      
      ;; Write out bias-subtracted image
      hires[indx[q]].img_ov = outfil
      hires[indx[q]].flg_ov = 1

      ;; Header
      mkhdr, main_head, ovimg
      sxdelpar, main_head, 'END'
      sxdelpar, head0, 'NAXIS'
      sxdelpar, head0, 'BITPIX'
      sxdelpar, head0, 'BZERO'
      sxdelpar, head0, 'SIMPLE'
      sxdelpar, head0, 'DATSUM'
      sxdelpar, head0, 'DATE'

      ;; Single chip
      if strmatch(hires[0].ccd,'SINGLE') then begin
         sxdelpar, head0, 'NAXIS'
         sxdelpar, head0, 'NAXIS1'
         sxdelpar, head0, 'NAXIS2'
         sxdelpar, head0, 'BITPIX'
      endif

      nhd = n_elements(head0)
      for qq=0L,nhd-1 do begin
          if strlen(strtrim(head0[qq],2)) GT 0 then $
            main_head = [main_head, head0[qq]]
      endfor
      
      mwrfits, ovimg, outfil, main_head, /create, /silent
      
  endfor

  print, 'hires_subbias: All done!'

  return
end
