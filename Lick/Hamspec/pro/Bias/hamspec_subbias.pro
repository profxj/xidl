;+ 
; NAME:
; hamspec_subbias   
;     Version 1.2
;
; PURPOSE:
;    Remove bias from all exposures listed in the "hamspec" structure.
;    The main driver is hamspec_suboscan which strips the image of the
;    overscan region (see that program for a full
;    description).  The hamspec_subbias routine will also remove an
;    archived bias image if requested.
;
; CALLING SEQUENCE:
;   
;  hamspec_subbias, hamspec, indx, /usebias, /nobiasrow, /clobber, /ARC,
;  /debug, BADROW=, RBIN=, CBIN=
;
; INPUTS:
;   hamspec  -  MIKE structure
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
;             were used in hamspec_mkbias (i.e. /nobiasrow )
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
;   hamspec_subbias, hamspec, indx, /usebias, /nobiasrow
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
;   hamspec_subbias_sngl -- Subroutine under hamspec_subbias.  Uses most of
;                        the keywords described above.  Example:  
;                        rslt = hamspec_subbias_sngl('Raw/mb0020.fits')
;   hamspec_suboscan
;
; REVISION HISTORY:
;   01-Feb-2005 Written by JXP  (taken from mike_subbias)
;                                
;------------------------------------------------------------------------------

function hamspec_subbias_sngl, rawfil, NAMP=namp, $
                             USEBIAS=usebias, NOBIASROW = nobiasrow, $
                             CLOBBER=clobber, DEBUG=debug, IMTYP=imtyp, $
                             BADROW=badrow, OVIMG=ovimg, NOFITS=nofits, $
                             SILENT=silent, CBIN=cbin, RBIN=rbin, FRAME=frame, $
                               GAIN_RATIO=gain_ratio

  colors = GetColor(/Load, Start=1)

  if  N_params() LT 1  then begin 
      print,'Syntax:  ' + $
        'rslt = hamspec_subbias_sngl(rawfil, /noBIASROW, ' + $
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

  if not keyword_set(NAMP) then stop
  if not keyword_set( FRAME ) then begin
     head = xheadfits(rawfil)
     frame = sxpar(head, 'OBSNUM')
      ;if pos1 NE -1 then stop
      ;pos = strpos(rawfil, '.fits')
      ;frame = long(strmid(rawfil,pos-4,4))
  endif

  ;; Check for output
  outfil = hamspec_getfil('ov_fil', FRAME=frame, /name, CHKFIL=chkf)
  if CHKF NE 0 and not keyword_set( CLOBBER ) then begin
      print,$
        'hamspec_subbias_sngl: File ', outfil, $
        ' found. Overscan already subtracted.'
      return, 0
  endif
          
  ;; Open Raw image
  raw = xmrdfits(strtrim(rawfil,2),0, head, /silent, /fscale)
  head0 = head
  if not keyword_set( SILENT ) then $
    print, 'hamspec_subbias: Subtracting the overscan for image ', rawfil
  sz = size(raw, /dimensions)
      
  ;; Set cbin, rbin
  ccd = sxpar(head0,'DSENSOR')
  if strlen(ccd) EQ 0 then ccd = sxpar(head0,'CCDID')
  ccd = strcompress(ccd,/remove_all)

  case ccd of 
     'Loral2Kx2K' : begin
        if not keyword_set( CBIN ) then cbin = round(2048. / sz[0] )
        if not keyword_set( RBIN ) then rbin = round(4096. / sz[1] )
     end
     'e2vCCD203-824kx4kthin' : begin
        if not keyword_set( CBIN ) then cbin = round(4096. / sz[0] )
        if not keyword_set( RBIN ) then rbin = round(4136. / sz[1] )
     end
     else: stop
  endcase
      
  ;; JXP :: Am saving the entire science image now
  cover = sxpar(head, 'COVER')
  n1 = sxpar(head,'NAXIS1')
  n2 = sxpar(head,'NAXIS2')
  dats = [0,n1-cover*namp-1,0,n2-1]
  fincol = dats[1]
 
  ;; TRIM Vignetting
;  dats[3] = dats[3] < round(3990./rbin)

  ;; TRIM Red chip
;  if chip EQ 3L then dats[0] = dats[0] > round(110./cbin)
          
  if keyword_set (USEBIAS) then begin 
      bias = hamspec_getfil('bias_fil', SZ=[2048L/cbin,4096L/rbin],  FIL_NM=bias_fil)
      if not keyword_set( SILENT ) then $
        print, 'hamspec_subbias: Using BIAS file: ', bias_fil
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
         stop ;; Not setup for this
      ovimg = temporary(ovimg) - bias
      sxaddpar, head, 'BIAS', 'T', bias_fil
  endif

  ;; Ovesrcan
  ;x_suboscan, raw, head, ovimg, fincol+5, IMTYPE=imtyp, $
  ;  DEBUG= debug, SVBAD=badrow, CBIN=cbin, RBIN=rbin
  ;; Subtract overscan
  case NAMP of
     1: begin ;; Single amp
        x_suboscan, raw, head, ovimg, fincol+5, IMTYPE=imtyp, $
                    DEBUG= debug, SVBAD=badrow, CBIN=cbin, $
                    RBIN=rbin
        ovimg = ovimg[dats[0]:dats[1],dats[2]:dats[3]]
     end
     2: begin  ;; two amp
        ovimg = fltarr(dats[1]+1, dats[3]+1)
        sz = size(ovimg, /dimen)
        ;; Left
        x_suboscan, raw[0:dats[1]+cover-1,*], $
                    head, ovimgl, fincol+3, IMTYPE=imtyp, $
                    DEBUG= debug, SVBAD=badrow, CBIN=cbin, $
                    RBIN=rbin
        if not keyword_set(GAIN_RATIO) then gain_ratio = 1.
        ovimg[0:sz[0]/2-1,*] = ovimgl[0:sz[0]/2-1,*] / GAIN_RATIO
        
        ;; Right
        x_suboscan, raw, $
                    head, ovimgr, fincol+cover+1, IMTYPE=imtyp, $
                    DEBUG= debug, SVBAD=badrow, CBIN=cbin, $
                    RBIN=rbin
        ovimg[sz[0]/2:*,*] = ovimgr[sz[0]/2:sz[0]-1,*]
     end
     else: stop
  endcase
      
  ;; Trim
;  dats = dats-1
  ;ovimg = ovimg[dats[0]:dats[1],dats[2]:dats[3]]
  sxaddpar, head0, 'TRIM', 'T'
  sxaddpar, head0, 'ODSEC', strjoin(string(dats))
      
  ;; Rotate
  ovimg = transpose(ovimg)

  ;; Header
  mkhdr, main_head, ovimg
  sxdelpar, main_head, 'END'
  sxdelpar, head0, 'NAXIS'
  sxdelpar, head0, 'NAXIS1'
  sxdelpar, head0, 'NAXIS2'
  sxdelpar, head0, 'BITPIX'
  sxdelpar, head0, 'BZERO'
  sxdelpar, head0, 'SIMPLE'
  sxdelpar, head0, 'DATSUM'
  sxdelpar, head0, 'DATE'

  nhd = n_elements(head0)
  for qq=0L,nhd-1 do begin
      if strlen(strtrim(head0[qq],2)) GT 0 then $
        main_head = [main_head, head0[qq]]
  endfor
      
  ;; Write out bias-subtracted image
  if not keyword_set(NOFITS) then $
    mwrfits, ovimg, outfil, main_head, /create, /silent
      
  if not keyword_set( SILENT ) then $
    print, 'hamspec_subbias_sngl: All done!'

  return, 1
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hamspec_subbias, hamspec, indx, USEBIAS=usebias, NOBIASROW = nobiasrow, $
                  CLOBBER=clobber, DEBUG=debug, BADROW=badrow, GAIN_RATIO=gain_ratio
              ; , VIEW=view

  colors = GetColor(/Load, Start=1)

  if  N_params() LT 2  then begin 
      print,'Syntax:  ' + $
        'hamspec_subbias, hamspec, indx, /noBIASROW, /USEBIAS, OVROOT=, ' + $
        '/CLOBBER, /DEBUG, BADROW= [v1.1]'
;      print,'Recommended:  ' + $
;        'hamspec_subbias, hamspec, /noBIASROW, /USEBIAS, /CLOBBER '
      return
  endif 
  
  
  cbin = 0L
  rbin = 0L
  
  nindx = n_elements(indx)
  
  ;; Loop
  for q=0,nindx-1 do begin
      ;; Check for output
      outfil = hamspec_getfil('ov_fil', $
                   FRAME=hamspec[indx[q]].frame, /name, CHKFIL=chkf)
      if CHKF NE 0 and not keyword_set( CLOBBER ) then begin
          print,$
            'hamspec_subbias: File ', outfil, $
            ' found. Overscan already subtracted.'
          hamspec[indx[q]].img_ov = outfil
          hamspec[indx[q]].flg_ov = 1
          continue
      endif
          
      ;; Open Raw image
      raw = xmrdfits(strtrim(hamspec[indx[q]].rootpth+hamspec[indx[q]].img_root,2), 0, head, $
                     /silent, /fscale)
      head0 = xheadfits(strtrim(hamspec[indx[q]].rootpth+hamspec[indx[q]].img_root,2))
      print, 'hamspec_subbias: Subtracting the overscan for image ', $
        hamspec[indx[q]].rootpth+hamspec[indx[q]].img_root,  ' ', hamspec[indx[q]].type

      ;; FINCOL
      cover = sxpar(head, 'COVER') 
      n1 = sxpar(head,'NAXIS1')
      n2 = sxpar(head,'NAXIS2')
      dats = [0,n1-cover*hamspec[indx[q]].amp-1,0,n2-1]
      fincol = dats[1]

      ;; TRIM Vignetting
;      dats[3] = dats[3] < round(3990./hamspec[indx[q]].rowbin)

      ;; TRIM Red chip
;      if chip EQ 3L then dats[0] = dats[0] > round(110./hamspec[indx[q]].colbin)

      if keyword_set (USEBIAS) then begin 
         stop ;; Not doing this so far
          bias = hamspec_getfil('bias_fil', SZ=[2048L/cbin,4096L/rbin], FIL_NM=bias_fil)
          print, 'hamspec_subbias: Using BIAS file: ', bias_fil
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
         stop ;; Not setup for this
          ovimg = temporary(ovimg) - bias
          hamspec[indx[q]].flg_anly = 3
          sxaddpar, head0, 'BIAS', 'T', bias_fil
      endif

      ;; Subtract overscan
      case hamspec[indx[q]].amp of
         1: begin ;; Single amp
            x_suboscan, raw, head, ovimg, fincol+5, IMTYPE=hamspec[indx[q]].type, $ 
                        DEBUG= debug, SVBAD=badrow, CBIN=hamspec[indx[q]].colbin, $
                        RBIN=hamspec[indx[q]].rowbin
            ovimg = ovimg[dats[0]:dats[1],dats[2]:dats[3]]
         end
         2: begin  ;; two amp
            ovimg = fltarr(dats[1]+1, dats[3]+1)
            sz = size(ovimg, /dimen)
            ;; Left
            x_suboscan, raw[0:dats[1]+cover-1,*], $
                        head, ovimgl, fincol+3, IMTYPE=hamspec[indx[q]].type, $ 
                        DEBUG= debug, SVBAD=badrow, CBIN=hamspec[indx[q]].colbin, $
                        RBIN=hamspec[indx[q]].rowbin
            ;; Left
            x_suboscan, raw, $
                        head, ovimgr, fincol+cover+1, IMTYPE=hamspec[indx[q]].type, $ 
                        DEBUG= debug, SVBAD=badrow, CBIN=hamspec[indx[q]].colbin, $
                        RBIN=hamspec[indx[q]].rowbin
            
            ;; Save
            if not keyword_set(GAIN_RATIO) then begin
               ;;gain_kludge = 1.0732  ;; JXP 19 feb 2013
               rtio = ovimgl[sz[0]/2-1,*] / ovimgr[sz[0]/2,*]
               gain_ratio = median(rtio)
               print, 'hamspec_subbias: Gain ratio measured at ', gain_ratio
            endif 

            ovimg[0:sz[0]/2-1,*] = ovimgl[0:sz[0]/2-1,*] / GAIN_RATIO
            ovimg[sz[0]/2:*,*] = ovimgr[sz[0]/2:sz[0]-1,*]
         end
         else: stop
      endcase
      
      ;; Trim
;      dats = dats-1
      head0 = head
      sxaddpar, head0, 'TRIM', 'T'
      sxaddpar, head0, 'ODSEC', strjoin(string(dats))

      ;; Rotate
      ovimg = transpose(ovimg)
      
      ;; Write out bias-subtracted image
      hamspec[indx[q]].img_ov = outfil
      hamspec[indx[q]].flg_ov = 1

      ;; Header
      mkhdr, main_head, ovimg
      sxdelpar, main_head, 'END'
      sxdelpar, head0, 'NAXIS'
      sxdelpar, head0, 'NAXIS1'
      sxdelpar, head0, 'NAXIS2'
      sxdelpar, head0, 'BITPIX'
      sxdelpar, head0, 'BZERO'
      sxdelpar, head0, 'SIMPLE'
      sxdelpar, head0, 'DATSUM'
      sxdelpar, head0, 'DATE'

      nhd = n_elements(head0)
      for qq=0L,nhd-1 do begin
          if strlen(strtrim(head0[qq],2)) GT 0 then $
            main_head = [main_head, head0[qq]]
       endfor
      
      mwrfits, ovimg, outfil, main_head, /create, /silent
      
  endfor

  print, 'hamspec_subbias: All done!'

  return
end
