;+ 
; NAME:
; ism_msngcog
;  V1.1
;
; PURPOSE:
;  Given a set of wavelengths, ew and sigew, calculate N,b for
;  a set of ions
; CALLING SEQUENCE:
;   
;  ism_msngcog, all_cog, NLMT=Nlmt, BLMT=blmt, /CHICHK, PLTONLY=
;    NSTP=, BSTP=, PSFILE=, OUTFIL=, /EXACT, ZLBL=
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   OUTFIL -- File with best fit values and error
;   PSFILE -- File for postscript plot
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   15-Apr-2006 Written by JXP  (based on ism_sngcog)
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro ism_msngcog, all_cog, NLMT=Nlmt, BLMT=blmt, CHICHK=chichk, $
                 NSTP=nstp, BSTP=bstp, PLOT=plot, STRCT=strct, $
                 DCHISQ=dchisq, NARR=narr, BARR=barr, SILENT=silent, $
                 GUESS_OFF=guess_off

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
      'ism_sngcog, wrest, ew, sigew, [Nlmt, blmt], /CHICHK, /PLOT,' + $
      'NSTP=, BSTP= (v1.1)' 
    return
  endif 

  if not keyword_set( GUESS_OFF ) then $
    guess_off = replicate(0., n_elements(all_cog))
  if not keyword_set( NLMT ) then nlmt = [10., 18.]
  if not keyword_set( BLMT ) then blmt = [5., 50.]
  if not keyword_set( NSTP ) then nstp = 100L
  if not keyword_set( BSTP ) then bstp = 100L
  if not keyword_set(EXACT) then begin
      print, 'fuse_cog: Using spline interpolation for the COG.  '
      print, '      This assumes a Maxwellian profile'
      cogmax_fil = getenv('XIDL_DIR')+'/Spec/Analysis/cogmax_tab.fits'
      cog_strct = xmrdfits(cogmax_fil, 1, /silent)
  endif

  nlvl = n_elements(all_cog)
  strct = { $
            N: fltarr(nlvl), $
            b: 0., $
            Nsig: fltarr(nlvl,2), $
            bsig: fltarr(nlvl,2) $
          }
  
  ;; CHISQ array
  all_chisq = dblarr(NSTP,BSTP,nlvl)
  tot_chi = dblarr(BSTP)

  ;; b_eff
  bval = blmt[0] + (blmt[1]-blmt[0])*findgen(bstp)/float(bstp)

  ;; LOOP
  for mm=0L,nlvl-1 do begin
      ;; Call sngcog
      ism_sngcog, all_cog[mm], NLMT=Nlmt+GUESS_OFF[mm], BLMT=blmt, $
                  CHISQ=chisq, PLOT=plot, NSTP=nstp, BSTP=bstp, /NOSIG
      
      ;; Save for Error Analysis
      all_chisq[*,*,mm] = chisq

      ;; Collapse and Add
      tot_chi = tot_chi + min(chisq, dimension=1)
  endfor

  if keyword_set( CHICHK ) then x_splot, bval, tot_chi, /bloc
      
  ;; Find min chisq
  mnchi = min(tot_chi, ichi)
  bgd = bval[ichi]
  strct.b = bgd

  ;; Get the best N values
  for mm=0L,nlvl-1 do begin
      ;; Nvalues
      Nval = Nlmt[0] + (Nlmt[1]-Nlmt[0])*findgen(nstp)/float(nstp) + GUESS_OFF[mm]
      mn = min(all_chisq[*,ichi,mm], imn)
      strct.N[mm] = Nval[imn]
  endfor
  return
      
  ii = ichi MOD NSTP
  jj = (ichi/NSTP)
  Ngd = Nval[ii]

;  print, 'Nbest = ', Ngd, ' bbest = ', bgd
  strct.N = Ngd
  strct.b = bgd

  ;; Plot
  dchisq = chisq - mnchi
  if keyword_set(PLOT) then begin
      clr = getcolor(/load)
      contour, dchisq, Narr, barr, levels=[2.3, 6.17, 18.4], $
            c_annota=['1!9s!X', '2!9s!X', '3!9s!X'], color=clr.black, $
            background=clr.white, chars=1.5, c_charsi=1.5
  endif
  

  ;; Calculate the probability array
  prob_arr = 1. - chisqr_pdf(chisq, 2)
  tot_prob = total(prob_arr)
  
  ;; Error in N
  sz = size(prob_arr)
  for i=0L, 2*NSTP/3 do begin
      if ii-i ge 0 and ii+i lt sz[1] then begin
          sub_prob = total( prob_arr[(ii-i)>0:ii+i,*] )
          if sub_prob GT 0.683*tot_prob then break
      endif else begin
          print,'ism_sngcog: cannot determine sigN'
          sigN = -9.99
          goto,assume_sigN
      endelse 
  endfor
  if i EQ 0 then stop
  sigN = (Nlmt[1]-Nlmt[0])*float(i)/float(nstp)
  assume_sigN:
  
  ;; Error in b
  for j=0L, 2*BSTP/3 do begin
      if jj-j ge 0 and jj+j lt sz[2] then begin
          sub_prob = total( prob_arr[*,(jj-j)>0:jj+j] )
          if sub_prob GT 0.683*tot_prob then break
      endif else begin
          print,'ism_sngcog: cannot determine sigb'
          sigb = -9.99
          goto,assume_sigb
      endelse
  endfor
  if j EQ 0 then stop
  sigb = (blmt[1]-blmt[0])*float(j)/float(bstp)
  assume_sigb:
  
  if not keyword_set(SILENT) then begin
      print, 'ism_sngcog: N = ', Ngd, '+/-', sigN
      print, 'ism_sngcog: b = ', bgd, '+/-', sigb
  endif

  return
end
