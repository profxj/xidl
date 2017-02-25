;+ 
; NAME:
; imacsls_skysub   
;     Version 1.1
;
; PURPOSE:
;    Sky Subtract image and add sky subtracted frame to the final
;    image (e.g. Final/f_imacsls0020.fits).  The program does a simply
;    row by row sky subtraction using x_fitrej.
;
; CALLING SEQUENCE:
;  imacsls_skysub, imacsls, setup, obj_id, side, [exp], 
;     /CHK, /STD, IAPER=, ORDR=, /NOVAC
;
; INPUTS:
;   imacsls -  IMACS long slit structure
;   setup   -  Setup ID value
;   obj_id  -  Object ID  (e.g. 0L, 1L, etc)
;   [exp]   -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /STD     - Sky subtract standard star
;  /CHK     - Show steps along the way 
;  ORDR=    - Order of sky POLY fit [default: 1]
;  /CLOBBER - Overwrite any previos sky image
;  IAPER    - Size of input aperture [default: use .aper tag]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   imacsls_skysub, imacsls, 1L, [0L], /CHK, ORDR=7L   
;           {Sky sub exposure 0 and order 7 only}
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_fitrej
;
; REVISION HISTORY:
;   09-Dec-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro imacsls_skysub, imacsls, setup, obj_id, side, exp, CHK=chk, STD=std, $
                    IAPER=iaper, $
                    USEOLD=useold, SEDG_FIL=sedg_fil, $
                    BCHK=bchk, indx=indx, X1=x1, x2=x2

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'imacsls_skysub, imacsls, setup, obj_id, side, [exp], /CHK, /STD, , '
      print, '      ORDR=, /USEOLD, /BCHK [v1.0]'
      return
  endif 
  
;  Optional Keywords

;  Find all relevant obj
  if not keyword_set( STD ) then begin
      indx = where(imacsls.flg_anly NE 0 AND imacsls.side EQ side AND $
                   imacsls.obj_id EQ obj_id AND imacsls.setup EQ setup AND $
                   strtrim(imacsls.type,2) EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'imacsls_skysub: No images to sky subtract!', obj_id
          return
      endif
  endif else begin
      indx = obj_id[0]
      nindx = 1L
  endelse
      
;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

  if not keyword_set( x1 ) then begin
      case side of
          1: begin
              x1 = round(1118. / imacsls[indx[0]].cbin)
              x2 = round(1235. / imacsls[indx[0]].cbin)
          end
          2: begin
              x1 = round(810.*2. / imacsls[indx[0]].cbin)
              x2 = round(924.*2. / imacsls[indx[0]].cbin)
          end
          else: stop
      endcase
  endif

;  Loop

  for q=0L,n_elements(exp)-1 do begin

      ;; Open Obj file
      objfil = imacsls_getfil('obj_fil', subfil= imacsls[indx[exp[q]]].img_root,$
                              /name, CHKF=chkf)
      if CHKF EQ 0 then begin
          print, 'imacsls_skysub: No Obj file! ', objfil, ' Skipping...'
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)

      ;; Aper
      if not keyword_set( IAPER ) then aper=objstr.aper else aper=iaper

      ;; IMG+VAR Fil 
      imgfil = strtrim(imacsls[indx[exp[q]]].img_final,2)
      if x_chkfil(imgfil+'*') EQ 0 then begin
          print, 'imacsls_skysub: Image file doesnt exist!', imgfil
          return
      endif
      print, 'imacsls_skysub: Reading Image files... ', imgfil
      img = xmrdfits(imgfil, 0, head, /silent)
      ivar = xmrdfits(imgfil, 1, /silent)
      sz_img = size(img, /dimensions)

      img_new = fltarr(sz_img[0],sz_img[1])

      sky = fltarr(x2-x1+1,sz_img[1])
      ;; Loop on rows
      nreg = 2L

      ; Set Regions
      nwreg = fltarr(2,2)
      nwreg[0,0] = 0
      nwreg[1,1] = (x2-x1)
      ;; Fit
      skyfstr = { fitstrct }
      skyfstr.func = 'POLY'
      if not keyword_set( ORDR ) then skyfstr.nord = 1 else skyfstr.nord=ordr
      skyfstr.lsig = 2.
      skyfstr.hsig = 2.
      skyfstr.niter = 5
      skyfstr.flg_rej = 1
      skyfstr.maxrej = 10
      if keyword_set( STD ) then skyfstr.nord = 0

      print, 'imacsls_skysub: Fitting...'
      dumx = findgen((x2-x1)+1)
      for jj=0L,sz_img[1]-1 do begin
          if keyword_set( STD ) then begin
              nwreg[0,0] = objstr.trace[jj]-12.-x1
              nwreg[0,1] = objstr.trace[jj]-aper[0]-x1
              nwreg[1,0] = objstr.trace[jj]+aper[1]-x1
              nwreg[1,1] = objstr.trace[jj]+12-x1
          endif else begin
              nwreg[0,1] = objstr.trace[jj]-aper[0]-x1
              nwreg[1,0] = objstr.trace[jj]+aper[1]-x1
          endelse
          sky[*,jj] = x_fitrej(dumx,img[x1:x2,jj], FITSTR=skyfstr, REG=nwreg,$
                                   IVAR=ivar[x1:x2,jj]) 
      endfor
      img_new[x1:x2,*] = img[x1:x2,*] - sky

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Ouptut New Image
      print, 'imacsls_skysub: Writing output to: ', imgfil
      if keyword_set( CHK ) then xatv, img_new, /block, min=-20, max=100.
      mwrfits, img, imgfil, head, /create, /silent
      mwrfits, ivar, imgfil, /silent
      mwrfits, img_new, imgfil, /silent
      ;; COMPRESS
      print, 'imacsls_skysub: Compressing...'
      spawn, 'gzip -f '+imgfil
  endfor
  
;  DONE
  print, 'imacsls_skysub: All done! '
  return
end
