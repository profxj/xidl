;+ 
; NAME:
; esi_subbias   
;     Version 1.1
;
; PURPOSE:
;    Median combine all ZRO frames (darks)
;      WARNING!  Assumes images are all of 1 mode (e.g. IMG, ECH, LWD)!!
;
; CALLING SEQUENCE:
;   
;  esi_subbias, esi, indx
;
; INPUTS:
;   esi   -  ESI structure
;   indx  -  Index numbers of frame to subtract (default output is OV)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  BIASFIL= - Name of bias file (default: Bias/BiasS[I].fits)
;  OVROOT=  - Root name of OV file (default: OV/ov_ )
;  /FORCE   - Overwrite existing OV files 
;  /NOOSCAN - Raw file has no overscan region (e.g. MIRA failure)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;  Currently only good for 1x1 binning
;
; EXAMPLES:
;   esi_subbias, esi, [47L,48L,49L]
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Jul-2002 Written by JXP
;   01-Feb-2003 Polished (JXP)
;   17-May-2013 Modified to work with bad windowing in 2x1 MR
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_subbias, esi, indx, OVROOT=OVROOT, BIASFIL=biasfil, FORCE=force $
                 , NOBAD = nobad, NOHOT = NOHOT, HOT_THRESH = HOT_THRESH $
                 , NOOSCAN=nooscan

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_subbias, esi, indx, OVROOT=, BIASFIL=, /FORCE [v1.1]'
      return
  endif 
  rbin = esi[indx[0]].rbin
  cbin = esi[indx[0]].cbin

  
;  Optional Keywords
  if not keyword_set( OVROOT ) then ovroot = 'OV/ov_'
  if not keyword_set( BIASFIL ) then begin
      biasfil = esi_getfil('bias_fil', esi[indx[0]].mode, $
                           cbin=cbin, rbin=rbin, /name)
  endif

  
; Open Bias file
  print, 'esi_subbias: Using BIAS file: ', biasfil
  if x_chkfil(biasfil+'*') EQ 0 then begin
      print, 'esi_subbias: Bias file not found.  Create first!', biasfil
      stop
  endif
  bias = xmrdfits(biasfil, /silent)
  sz = size(bias, /dimensions)

; Loop

  for q=0,n_elements(indx)-1 do begin
      ;; Check binning
      if esi[indx[q]].rbin NE rbin OR esi[indx[q]].cbin NE cbin then begin
          print, 'esi_subbias:  Binning varies in these files!  Trouble..'
          stop
      endif
      ;; Check for output
      outfil = ovroot+esi[indx[q]].img_root
      a = findfile(outfil+'*', count=na)
      if na NE 0 and not keyword_set( FORCE ) then begin
          print, 'esi_subbias: File ', outfil, ' found. Not resubtracting'
          esi[indx[q]].img_ov = outfil
          esi[indx[q]].flg_ov = 1
          continue
      endif
      ;; Open Raw image
      raw = xmrdfits(esi[indx[q]].rootpth+esi[indx[q]].img_root, 0, head, $
                    /silent, /fscale)
      print, 'esi_subbias: Subtracting the bias for ', $
        esi[indx[q]].rootpth+esi[indx[q]].img_root

      ;; Check for bad windowing - occurs if windowing reset during
      ;; observations if ESI has issues - MR 2013
      ;; Also deal if there is *NO* overscan on science image (but ok BIAS) -- JXP 2014

      rawsize = size(raw)
      badwindow = 0
      if rawsize[1] eq 2208 then badwindow = 1
      if rawsize[1] eq 2048 then badwindow = 2 

      ;; Allow for mode
      case esi[indx[q]].mode of
          0: begin
             if keyword_set(NOOSCAN) then stop
              if esi[indx[q]].cbin NE 1 OR esi[indx[q]].rbin NE 1 then stop
              ovimg = raw[160:910,80:1425L]  - bias ; IMG
          end
          1: begin ; LOWD
             if keyword_set(NOOSCAN) then stop
              if esi[indx[q]].cbin NE 1 OR esi[indx[q]].rbin NE 1 then stop
              if sz[0] LT 2100L then stop
              case esi[indx[q]].namp of
;                  1: ovimg = raw[12:2059,230:3529L] - bias[*,230:3529L] 
;                  2: ovimg = raw[25:2072,230:3529L] - bias[*,230:3529L]
                  1: ovimg = raw[25:2072,230:3529L] - bias[25:2072,230:3529L] 
                  2: ovimg = raw[25:2072,230:3529L] - bias[25:2072,230:3529L]
                  else: stop
              endcase
          end
          2: begin ; ECH
             ovimg = fltarr(2048L/esi[indx[q]].cbin,sz[1])
;              if sz[0] LT 2100L then stop
             if keyword_set(NOOSCAN) then begin
                ;; Occasional trimming of image by MIRA algorithm
                ovimg = raw - bias[0:(2048L/esi[indx[q]].cbin)-1,*]
             endif else begin
                ;; ???? This badwindow stuff is all hard-wired for 1x1 binning
                ;; bad window case MR 2013
                if badwindow EQ 1 THEN begin
                   tmp = raw - bias[24L:*, *]
                endif else if badwindow EQ 2 THEN begin
                   tmp = raw - bias[24L:(2048L+24L-1L), *]
                endif else begin
                   tmp = raw - bias
                endelse
                case esi[indx[q]].namp of
;                  1: ovimg = raw[12:2059,*] - bias 
;                  2: ovimg = raw[25:2072,*] - bias 
                   1: begin
                      if cbin NE 1 then stop
                      ;; Full chip
                      ov1 = djs_median(tmp[2072:2133,*],1)
                      fitstr = x_setfitstrct()
                      f1 = x_fitrej(findgen(4096L/rbin), ov1,fitstr=fitstr)
                      ovimg[*,*] = tmp[12:2059,*] - replicate(1.,2048L)#f1
                   end
                   2: begin 
                      case cbin of 
                         1: begin
                            ;; Bad window case - Lose overscan on left
                            ;; which is 24 pixels wide. MR 2013
                            if badwindow EQ 1 THEN BEGIN 
                               fixwin = 24L
                               x11 = 25-fixwin
                               x12 = 1047L-fixwin
                               x21 = 1048L-fixwin
                               x22 = 2072L-fixwin
                               a11 = 2083L-fixwin
                               a12 = 2137L-fixwin
                               a21 = 2160L-fixwin
                               a22 = 2220L-fixwin
                            endif else begin
                               x11 = 25
                               x12 = 1047L
                               x21 = 1048L
                               x22 = 2072L
                               a11 = 2083L
                               a12 = 2137L
                               a21 = 2160L
                               a22 = 2220L
                            endelse
                         end
                         2: begin
                            x11 = 12
                            x12 = 523L
                            x21 = 524L
                            x22 = 1035L
                            a11 = 1050L
                            a12 = 1072L
                            a21 = 1078L
                            a22 = 1105L
                         end
                         4: begin
                            x11 = 6
                            x12 = 261L
                            x21 = 262L
                            x22 = 517L
                            a11 = 520L
                            a12 = 536L
                            a21 = 539L
                            a22 = 556L
                         end
                         else: stop
                      endcase
                      ;; For frames for which there is no overscan on
                      ;; science frames but an OK bias, we simply
                      ;; subtracted the bias above. No overscan
                      ;; subtraction JFH 03.2016
                      IF badwindow EQ 2 THEN ovimg = tmp $
                      ELSE BEGIN 
                         ;; LHS
                         ov1 = djs_median(tmp[a11:a12, *], 1)
                         fitstr = x_setfitstrct()
                         f1 = x_fitrej(findgen(4096L/rbin), ov1, fitstr = fitstr)
                         ovimg[0:x12-x11, *] = $
                            tmp[x11:x12, *] - replicate(1., x12-x11+1)#f1
                         ;; RHS
                         ov2 = djs_median(tmp[a21:a22, *], 1)
                         fitstr = x_setfitstrct()
                         f2 = x_fitrej(findgen(4096L/rbin), ov2, fitstr = fitstr)
                         ovimg[x12-x11+1:*, *] = tmp[x21:x22, *] - $
                                                 replicate(1., x22-x21+1)#f2
                      ENDELSE
                   END
                else: stop
             endcase
          endelse
          ;; Mask bad pix
          if not keyword_set( NOBAD ) then $
             esi_echbadpix, esi[indx[q]], ovimg, NOHOT = NOHOT $
                            , HOT_THRESH = HOT_THRESH
       end
      else: begin
         print, 'esi_subbias: Not prepared for mode ', esi[indx[q]].mode
         stop
      end
   endcase
  
  ;; Update flg and header
  esi[indx[q]].flg_anly = 3
  sxaddpar, head, 'BIAS', 'T', biasfil
  
  ;; Write
  esi[indx[q]].img_ov = outfil
  esi[indx[q]].flg_ov = 1
      mwrfits, ovimg, outfil, head, /create, /silent
   endfor

  print, 'esi_subbias: All done!'

  return
end
