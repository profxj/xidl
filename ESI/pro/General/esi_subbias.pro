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
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_subbias, esi, indx, OVROOT=OVROOT, BIASFIL=biasfil, FORCE=force, $
                 NOBAD=nobad

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_subbias, esi, indx, OVROOT=, BIASFIL=, /FORCE [v1.1]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( OVROOT ) then ovroot = 'OV/ov_'
  if not keyword_set( BIASFIL ) then begin
      if esi[indx[0]].mode EQ 0 then biasfil = 'Bias/BiasI.fits' $
        else biasfil = 'Bias/BiasS.fits'
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

      ;; Allow for mode
      case esi[indx[q]].mode of
          0: begin
              if esi[indx[q]].cbin NE 1 OR esi[indx[q]].rbin NE 1 then stop
              ovimg = raw[160:910,80:1425L]  - bias ; IMG
          end
          1: begin ; LOWD
              if esi[indx[q]].cbin NE 1 OR esi[indx[q]].rbin NE 1 then stop
              if sz[0] LT 2100L then stop
              case esi[indx[q]].namp of
;                  1: ovimg = raw[12:2059,230:3529L] - bias[*,230:3529L] 
;                  2: ovimg = raw[25:2072,230:3529L] - bias[*,230:3529L]
                  1: ovimg = raw[*,230:3529L] - bias[*,230:3529L] 
                  2: ovimg = raw[*,230:3529L] - bias[*,230:3529L]
                  else: stop
              endcase
          end
          2: begin  ; ECH
;              if sz[0] LT 2100L then stop
              tmp = raw - bias
;              ovimg = fltarr(2072L-25+1,4096L)
              ovimg = fltarr(2048L/esi[indx[q]].cbin,sz[1])
              case esi[indx[q]].namp of
;                  1: ovimg = raw[12:2059,*] - bias 
;                  2: ovimg = raw[25:2072,*] - bias 
                  2: begin
                      case esi[indx[q]].cbin of 
                          1: begin
                              x11 = 25
                              x12 = 1047L
                              x21 = 1048L
                              x22 = 2072L
                              a11 = 2083L
                              a12 = 2137L
                              a21 = 2160L
                              a22 = 2220L
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
                      ;; LHS
                      ov1 = djs_median(tmp[a11:a12,*],1)
                      fitstr = x_setfitstrct()
                      f1 = x_fitrej(findgen(4096L/esi[indx[q]].rbin), $
                                    ov1,fitstr=fitstr)
                      ovimg[0:x12-x11,*] = $
                        tmp[x11:x12,*] - replicate(1.,x12-x11+1)#f1
                      ;; RHS
                      ov2 = djs_median(tmp[a21:a22,*],1)
                      fitstr = x_setfitstrct()
                      f2 = x_fitrej(findgen(4096L/esi[indx[q]].rbin), $
                                    ov2, fitstr=fitstr)
                      ovimg[x12-x11+1:*,*] = tmp[x21:x22,*] - $
                        replicate(1.,x22-x21+1)#f2
                  end
                  else: stop
              endcase
              ;; Mask bad pix
              if not keyword_set( NOBAD ) then $
                esi_echbadpix, esi[indx[q]], ovimg
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
