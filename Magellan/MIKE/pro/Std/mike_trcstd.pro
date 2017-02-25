;+ 
; NAME:
; mike_trcstd   
;     Version 1.0
;
; PURPOSE:
;    Trace the std through each order of the image.  The program
;    uses a standard star as a crutch (/STD) through regions where the object
;    has very low flux.  The trace is written into the object
;    structure (e.g. Extract/Obj_mike0024.fits)
;
; CALLING SEQUENCE:
;   
;  mike_trcstd, mike, setup, side
;
; INPUTS:
;   mike     -  ESI structure
;   obj_id  -  Object ID  (e.g. 0L, 1L, etc)
;   [exp_id]   -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /CHK    - Show final trace
;   /STD    - Use standard star as a crutch
;   /FAINT  - Faint object; sum up more rows (40) to search for flux 
;   NCOLL=  - Set number of rows to sum by hand (default: 25)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_trcstd, mike, 1L, [0L], /CHK, /STD
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   22-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_trcstd, mike, setup, side, CHK=chk, $
                   SEDG_FIL=sedg_fil, DEBUG=debug, ORDRS=ordrs, FRAD=frad

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_trcstd, mike, setup, side,  ' + $
        '/CHK, NCOLL=, /GUIDE, /DEBUG [v1.0]'
      return
  endif 

;  Optional Keywords
  if not keyword_set( SIDE ) then side = [1L,2L]
  if not keyword_set( FRAD ) then frad = 10.
;  if keyword_set( USESTD ) or keyword_set( GUIDE ) then flg_guide = 1
  flg_guide = 1

  resolve_routine, 'mike_trcobj', /no_recompile

;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOP ON SIDE

  for ii=0L,n_elements(side)-1 do begin
      
      qq = side[ii]
      ;; INDX
      indx = where(mike.flg_anly NE 0 AND mike.side EQ qq AND $
                   mike.setup EQ setup AND $
                   strtrim(mike.type,2) EQ 'STD', nindx) 

      ;;  Setup and setup
      ordr_str = mike_getfil('ordr_str', setup, SIDE=qq)
      nordr = n_elements(ordr_str)

      ordrs = [ordr_str[0].order, ordr_str[nordr-1].order]

;;;;;;;;;;;;;;;;;;;;;;;;;
; Loop on indx

      for nn=0L,nindx-1 do begin
          idx = indx[nn]
          xyoff = mike[idx].arc_xyoff[0]

          ;; Edges
          ordr_str.lhedg = ordr_str.lhedg + xyoff
          ordr_str.rhedg = ordr_str.rhedg + xyoff

          ;; Create obj fil
          tmp = { dblsobjstrct }
          objstr = replicate(tmp, nordr)
          objstr.slit_fil = ' '
          objstr.spec2d_fil = ' '
          objstr.img_fil = ' '
          objstr.UT = ' '
          objstr.instr_strct = ' '
          objstr.field = ' '
          objstr.exp = mike[idx].exp
          objstr.obj_id = 'a'
          nobj = 0L

          ;; Open Image, Variance
          print, 'mike_trcstd: Opening image...'
          imgfil = strtrim(mike[idx].img_final,2)
          img = xmrdfits(imgfil, 0, head, /silent)
          ivar = xmrdfits(imgfil, 1, /silent)
          sz_img = size(img, /dimensions)

;;;;;;;;;;;;;;;;;;;;;
; LOOP ON ORDERS

          if qq EQ 1 then ostep = 1 else ostep = -1
          for kk=ordrs[0],ordrs[1],ostep do begin
              ;; Edge
              flg_edge = 0 
              ;; Index
              mm = where(ordr_str.order EQ kk, nmm)
              if nmm EQ 0 then stop else mm = mm[0]
              print, 'mike_trcstd: Tracing order ', string(kk,2)

              rlhe = round(ordr_str[mm].lhedg)
              rrhe = round(ordr_str[mm].rhedg) 

              ;; Find at center if the entire order is on, otherwise
                  ;; go up
              if rrhe[sz_img[1]/2] LT sz_img[0] then medj = sz_img[1]/2 $
              else begin
                  b = findgen(sz_img[1])
                  a = where(rrhe LT sz_img[0] AND b GT sz_img[1]/2, na)
                  if na EQ 0 then stop else medj=a[0]
              endelse

              ;; Find obj
              smsh = djs_median(img[rlhe[medj]:rrhe[medj],medj-2:medj+2],2)
              center = x_centspln(findgen(n_elements(smsh)), smsh, /FORCE)
              gdcen = rlhe[medj] + center
              objstr[qq].ycen = gdcen
              objstr[qq].xcen = float(medj)

              ;; Crude guess at trace
              frac = (gdcen - ordr_str[mm].lhedg[medj]) / $
                (ordr_str[mm].rhedg[medj]-ordr_str[mm].lhedg[medj])
              crude_trc = frac*(ordr_str[mm].rhedg-ordr_str[mm].lhedg)+ $
                ordr_str[mm].lhedg

              ;; TRACE
              if qq EQ 1 then ncoll = 5L else ncoll = 3L
              mike_trcobj_sngl, img, ivar, setup, qq, mm, trace, NCOLL=ncoll, $
                GUIDE=crude_trc, XYOFF=xyoff, FRAD=frad, DEBUG=debug
;              mike_trcobj_sngl, img, ivar, setup, side, $
;                qq, trace, NCOLL=ncoll, GUIDE=crude_trc, XYOFF=xyoff, $
;                FRAD=frad, DEBUG=debug, MXORDR=7L, /SUBSKY

              ;; Fill up objstr
              objstr[mm].trace[0:sz_img[1]-1] = trace

          endfor

          ;;  Set 3" aperture for standard stars, +-1.5" round up.

          slitwidth = median(ordr_str.rhedg - ordr_str.lhedg)
          aperhalf = long(1.5*slitwidth/5.0)+1
          objstr.aper = [aperhalf,aperhalf]
         
          
          ;; CHK
          if keyword_set( CHK ) then begin
              tmp = img
              for kk=ordrs[0],ordrs[1],ostep do begin
                  mm = where(ordr_str.order EQ kk, nmm)
                  if nmm EQ 0 then stop else mm = mm[0]
                  trc_msk = round(objstr[mm].trace[0:sz_img[1]-1]) $
                    + lindgen(sz_img[1])*sz_img[0]
                  tmp[trc_msk] = -1000
              endfor
              xatv, tmp, /block, min=-20, max=1000
          endif
          
          ;; OUTPUT
          objfil = mike_getfil('obj_fil', subfil=mike[idx].img_root, /name)
          print, 'mike_trcstd: Creating obj file ', objfil
          mwrfits, objstr, objfil, /create, /silent
          spawn, 'gzip -f '+objfil
      endfor
  endfor
  
  print, 'mike_trcstd: All done!'

  return
end
              
;              trc_fit = x_setfitstrct(NITER=4L, NORD=10L, FLGREJ=1L, HSIG=2.2, $
;                                      LSIG=2.2, FUNC='POLY')
;              if qq EQ 1 AND kk GE 106 then begin
;                  trc_fit.hsig = 3.0
;                  trc_fit.lsig = 3.0
;              endif
;              ;; Red side
;              if qq EQ 2 then begin
;                  trc_fit.hsig = 2.5
;                  trc_fit.lsig = 2.5
;                  trc_fit.nord = 13
;                  trc_fit.niter = 5
;              endif
              
              ;; NORMAL FIT
;              new_trc = x_fitrej(yfit[0:gdfit-1], xfit[0:gdfit-1], $
;                                 SIG=xsig[0:gdfit-1], FITSTR=trc_fit)
;              print, 'mike_trcstd: RMS = ', trc_fit.rms

              ;; Test dx
;              x_splot, yfit[0:gdfit-1], $
;                xfit[0:gdfit-1]-ordr_cen[round(yfit[0:gdfit-1]),mm], /block, PSYM1=1
              ;; Test dx/Dx
;              x_splot, yfit[0:gdfit-1], $
;                (xfit[0:gdfit-1]-ordr_cen[round(yfit[0:gdfit-1]),mm]) / $
;                (ordr_str[mm].rhedg[round(yfit[0:gdfit-1])] $
;                 -ordr_str[mm].lhedg[round(yfit[0:gdfit-1])]), /block, PSYM1=1

              ;; Test frac
;              x_splot, yfit[0:gdfit-1], $
;                (xfit[0:gdfit-1]-ordr_str[mm].lhedg[round(yfit[0:gdfit-1])]) / $
;                (ordr_str[mm].rhedg[round(yfit[0:gdfit-1])] $
;                 -ordr_str[mm].lhedg[round(yfit[0:gdfit-1])]), /block, PSYM1=1

                  ;; objstr (Take trace to jend)
;              objstr[mm].trace[0:jend] = x_calcfit(findgen(jend+1), FITSTR=trc_fit)
