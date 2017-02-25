;+ 
; NAME:
; xoverscan
;   Version 1.2
;
; PURPOSE:
;    Subtracts an overscan region from an image
;
; CALLING SEQUENCE:
;   
;   xoverscan, img, ccd, [medov], OVSEC=, OVROOT=, /INTER, /ERASE
;    OVIMG=, /NOFITS, /SILENT, /MEDONLY, /COMPRESS, TRIMSEC=
;
; INPUTS:
;   img        - Raw fits file
;   ccd        - CCD name (Tek5, SITe1, LRISR, LRISb, GEN, WFTek5,
;                SITe3)
;
; RETURNS:
;
; OUTPUTS:
;   nwimg      - Writes overscan subtracted image to OV directory
;                  with an ov_ prefix (default)
;
; OPTIONAL KEYWORDS:
;   func       - String defining function to fit with (POLY is
;                default)
;   ordr       - Order of the function
;   ovsec      - String defining overscan section
;   ovroot     - String for OV root ('OV/' is the default)
;   erase      - Erase pre-existing ov files
;   NOFITS     - Suppress output to fits file
;   MEDONLY    - Useful in updating distruct
;   /COMPRESS  - Spawn gzip on the output image
;   TRIMSEC    - Image region to trim to
;   
;
; OPTIONAL OUTPUTS:
;   medov      - Median of each OV image
;   OVIMG      - Float array of OV image
;
; COMMENTS:
;
; EXAMPLES:
;   xoverscan, 'ccd001.fits', 'Tek5'
;
;
; PROCEDURES CALLED:
;  xregtovec
;  djs_median
;  x1dfit
;
; REVISION HISTORY:
;   20-June-2001 Written by JXP
;   05-Dec-2001  Added OVIMG option
;-
;------------------------------------------------------------------------------

pro xoverscan, img, ccd, medov, FUNC=func, ORDR=ordr, OVSEC=ovsec, $
               OVROOT=ovroot, INTER=inter, ERASE=erase, OVIMG=ovimg, $
               NOFITS=nofits, SILENT=silent, MEDONLY=medonly, $
               TRIMSEC=trimsec, COMPRESS=compress

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'xoverscan, image, ccd, [medov], FUNC=,ORDR=,OVSEC=,OVROOT=,'
    print, '        ERASE=,/INTER, OVIMG=, /NOFITS, /SILENT, MEDONLY=  [V1.2]'
    return
  endif 

; Allow img to be a strarr of fits files or a float array!

  type = size(img,/type)
  case type of 
      7: begin                  ; Array of image names
          nimg = n_elements(img) 
          flg_img = 0
      end
      else: begin
          nimg = 1
          flg_img = 1
          img_nm = 'tmp.fits'
          NOFITS = 1  ; Not allowing fits output for the moment
      end
  endcase
  

;  Optional Keywords
  if not keyword_set( FUNC ) then    func = 'POLY'
  if keyword_set( OVSEC ) then    flgovsec    = 1 else flgovsec = 0
  if not keyword_set( OVROOT ) then ovroot = 'OV/'
  if not keyword_set( ORDR ) then ordr = 5
  if arg_present( medov ) then medov = fltarr(nimg)  ; Median array

; Start looping

  for qq=0,nimg-1 do begin
      ; READ IN IMAGE
      case flg_img of 
          0: begin  ; FITS files
              img_nm = strtrim(img[qq],2)
              print, 'xoverscan: OV subtracting ', img_nm
              data = xmrdfits( img_nm, 0, header, /fscale, /silent)
              sz = SIZE( data )
          end
          1: begin  ; Float array
              print, 'xoverscan: OV subtracting '
              data = img
              sz = SIZE( data )
          end
      endcase

      ; Median Only!
      if keyword_set( MEDONLY ) then begin
          if not arg_present( medov ) then begin
              print, 'xoverscan: If you want the median, you need to set medov!'
              return
          endif
          case ccd of 
              'GEN' : medov[qq] = median(data) ; General CCD 
              'Tek5' : medov[qq] = median(data) ; Tek5 for Direct Imaging
              'WFTek5' : medov[qq] = median(data) ; Tek5+WFCCD
              'SITe1' : medov[qq] = median(data) ; SITe1
              'SITe3' : medov[qq] = median(data) ; SITe3
              'LRISb' : medov[qq] = median(data) ; LRIS blue
              'LRISR' : begin ; LRIS red
                  ncolm = 846
                  medov[qq] = median(data[0:ncolm-1,*])
              end
              else : print, ccd, ' is not a valid CCD'
          endcase
          continue
      endif

      ; DEAL with FITS output
      if not keyword_set( NOFITS ) then begin
          lastslsh = strpos(img_nm,'/',/reverse_search)
          if lastslsh NE -1 then nopth = strmid(img_nm, lastslsh+1) $
          else nopth = img_nm
          ovout = strjoin([ovroot,'ov_',nopth])
          a = findfile(ovout,count=count)
          if (count NE 0 AND not keyword_set(ERASE)) then begin
              print, ovout, ' exists: NOT overwriting!'
              ; Calculate median if arg is present
              if arg_present( medov ) then begin
                  a = xmrdfits(ovout, /silent, /fscale)
                  medov[qq] = median(a)
              endif
              continue
          endif
      endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; CCD specific stuff

      flg_twoamp = 0
      offset = 0
      flg_bias = 0
      bias_ord = 10
      case ccd of 
          'GEN' : begin  ; General CCD (takes last 50 columns)
              if(flgovsec EQ 1) then ova = xregtovec(ovsec,sz) $
              else ova = [sz[1]-51, sz[1]-1, 0, sz[2]-1] ;   Last 50 columns
              ;;; ncolm ;;;
              ncolm = ova[0]
              ;;; TRIM ;;;
              trimsec = [0, sz[2]-1]
          end
          'Tek5' : begin ; Tek5 for Direct Imaging
              if(flgovsec EQ 1) then ova = xregtovec(ovsec,sz) $
              else ova = [2050, 2147, 0, 2047] ;   Overscan region for Tek5 (Imaging)
              ;;; ncolm ;;;
              ncolm = 2048
              ;;; TRIM ;;;
              trimsec = [0, 2046]
              flg_bias = 1
              if flgovsec EQ 0 then ovhead = '[2050:2147,0:2047]' else ovhead=ovsec
          end
          'WFTek5' : begin
              if(flgovsec EQ 1) then ova = xregtovec(ovsec,sz) $
              else ova = [2050, 2079, 0, 2047] ;   Overscan region for WFTek5
              ova[1] = (sz[1]-1) < ova[1]
              ;;; ncolm ;;;
              ncolm = 2048
              ;;; TRIM ;;;
              trimsec = [0, 2046]
              if flgovsec EQ 0 then ovhead = '[2050:2079,0:2047]' else ovhead=ovsec
              flg_bias = 1
          end
          'SITe1' : begin 
              if(flgovsec EQ 1) then ova = xregtovec(ovsec,sz) $
              else ova = [2052, 2147, 0, 2047] ;   Overscan region for SITe1 (IMG)
              ;;; ncolm ;;;
              ncolm = 2048
              ;;; TRIM ;;;
              trimsec = ova
              if flgovsec EQ 0 then ovhead = '[2052:2147,0:2047]' else ovhead=ovsec
              flg_bias = 1
          end
          'SITe3' : begin 
              if(flgovsec EQ 1) then ova = xregtovec(ovsec,sz) $
              else ova = [2052, 2147L, 0, sz[2]-1] ; Region for SITe3 on Swope (IMG)
              ;;; ncolm ;;;
              ncolm = 2048L
              ;;; TRIM ;;;
              trimsec = [0, 3149L]
              if flgovsec EQ 0 then ovhead = '[2052:2147,*]' else ovhead=ovsec
              ;; BIAS
              flg_bias = 1
              bias_msk = replicate(1,2048L) 
              bias_msk[770:840] = 0
              bias_ord = 17
          end
          'LRISR' : begin
             ; determine whether we accidentally did a partial chip readout.
             window=sxpar(header, 'WINDOW')

             if window ne '0,0,0,2048,2048' $
                and window ne '0,0,550,2048,1000' then $
                   print, 'ERROR: not sure what portion of chip was read out.'

             CASE window OF
                '0,0,0,2048,2048' : begin
                   final = fltarr(1630,2048) ; Trim
                   ;;;;;  LHS ;;;;
                   ova = [2104, 2154, 0, 2047] ;   Overscan for LHS of LRIS
                   ncolm = 846
                   ;;;;;  RHS ;;;;
                   ovb = [2179, 2229, 0, 2047] ;   Overscan for RHS of LRIS
                   ncolm2 = 1630

                   ovsec  = '[2104:2154,0:2047] [2179:2229,0:2047]'
                   ovhead = '[2104:2154,0:2047] [2179:2229,0:2047]'
                   
                   offset = 220
                   flg_twoamp = 1
                end
                '0,0,550,2048,1000' : begin
                   final = fltarr(1630,1000) ; Trim
                   ;;;;;  LHS ;;;;
                   ova = [2104, 2154, 0, 999] ;   Overscan for LHS of LRIS
                   ncolm = 846
                   ;;;;;  RHS ;;;;
                   ovb = [2179, 2229, 0, 999] ;   Overscan for RHS of LRIS
                   ncolm2 = 1630

                   ovsec  = '[2104:2154,0:999] [2179:2229,0:999]'
                   ovhead = '[2104:2154,0:999] [2179:2229,0:999]'
                   
                   offset = 220
                   flg_twoamp = 1
                end
                ENDCASE
          end
          'LRISb' : begin 
              if(flgovsec EQ 1) then ova = xregtovec(ovsec,sz) $
              else ova = [2069, 2147, 0, sz[2]-1] ;   Overscan region for LRISb
              ;;; ncolm ;;;
              ncolm = 2069
              ;;; TRIM ;;;
              if not keyword_set( TRIMSEC ) then trimsec = [0, sz[2]-1]
              if flgovsec EQ 0 then ovhead = '[2069:2147,*]' else ovhead=ovsec
          end
          else : print, ccd, ' is not a valid CCD'
      endcase

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;; ; OV SUBTRACT ;;;;;;;;;;;

;;;;;;; ONE AMP ;;;;;;;;
      if flg_twoamp EQ 0 then begin 
          ; Take the median with djs_median
          ovmed = djs_median(data[ova[0]:ova[1],ova[2]:ova[3]],1)
          ;;; FIT ;;;
          ovspl = x1dfit(ovmed, func=func, $
                         nord=ordr, INTER=inter, LSIG=3., HSIG=3.)

          ;;; BIAS ;;;
          if flg_bias NE 0 then begin
              if not keyword_set( SILENT ) then $
                print, 'xoverscan: Subtracting the bias'
              biasrow = data[offset:offset+ncolm-1,sz[2]-1] $
                - replicate(ovspl[sz[2]-1],ncolm)
              ; Fit
              if not keyword_set( BIAS_MSK ) then $
                bias_msk = replicate(1, n_elements(biasrow))
              bias_fit = x_setfitstrct( NITER=4, HSIG=2., LSIG=2., FLGREJ=1, $
                                      FUNC='BSPLIN', NORD=bias_ord)
              bias = x1dfit(biasrow, fitstr=bias_fit, MSK=bias_msk)
              ;; Replaced masked region with unfit values!
              a = where(bias_msk EQ 0, na)
              if na NE 0 then bias[a] = biasrow[a]
              ;; Subtract
              tmp = replicate(1., trimsec[1]-trimsec[0]+1)
              final = data[offset:offset+ncolm-1, $
                           trimsec[0]:trimsec[1]] - bias#tmp
          endif

          ;;; SUBTRACT ;;;
          tmp = replicate(1, ncolm)

          if flg_bias EQ 0 then final = data[offset:offset+ncolm-1,$
                                             trimsec[0]:trimsec[1]] - tmp#ovspl $
          else final = temporary(final) - tmp#ovspl

          ;;; Median ;;;
          if arg_present( medov ) then medov[qq] = median(final)

;;;;;;; TWO AMP ;;;;;;;;;;;;;;;;;
      endif else begin  

          ;;; LHS ;;;
          ; MEDIAN
          ovmed = djs_median(data[ova[0]:ova[1],ova[2]:ova[3]],1)
          ; FIT 
          print, 'Fit LHS first'
          ovspl = x1dfit(findgen(sz[2]), ovmed, func=func, $
                         nord=ordr, INTER=inter, LSIG=3., HSIG=3.) 
          ; SUBTRACT
          tmp = replicate(1, ncolm)
          final[0:ncolm-1,*] = data[offset:offset+ncolm-1,*] - tmp#ovspl
;          for i=0L,ncolm-1 do final[i,*] = data[i+offset,*] - ovspl

          ;;; RHS ;;;
          ; MEDIAN
          ovmed = djs_median(data[ovb[0]:ovb[1],ovb[2]:ovb[3]],1)
          ; FIT 
          print, 'Fit RHS next'
          ovspl = x1dfit(findgen(sz[2]), ovmed, func=func, $
                         nord=ordr, INTER=inter, LSIG=3., HSIG=3.) 
          ; SUBTRACT 
          tmp = replicate(1, ncolm2-ncolm+1)
          final[ncolm:ncolm2-1,*] = $
            data[offset+ncolm:offset+ncolm2-1,*] - tmp#ovspl
;          for i=ncolm,ncolm2-1 do final[i,*] = data[i+offset,*] - ovspl
      
          ;;; Median ;;;
          ; LHS only
          if arg_present( medov ) then medov[qq] = median(final[0:ncolm-1,*])

      endelse

      ;;; Output ;;;
      if not keyword_set( NOFITS ) then begin
          sxaddpar, header, 'OVSCAN', ovhead
          mwrfits, final, ovout, header, /silent, /create
          if keyword_set(COMPRESS) then spawn, 'gzip -f '+ovout
      endif
      
      ;;; Output Image ;;;
      if arg_present( OVIMG ) then ovimg = final
      
      ;;; Free the memory ;;;
      delvarx, final

  endfor

end

