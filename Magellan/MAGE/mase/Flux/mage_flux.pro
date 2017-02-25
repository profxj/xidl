; + 
; NAME:
; mage_flux
; Version 0.1
;
; PURPOSE:
;  Uses a standard star fit calclulated by mage_fitstd to flux the
;  orders of an observation.
;
; CALLING SEQUENCE:
;
;  mage_flux,calsavfile,objstr,REJ=rej
;
; INPUTS:
;   calsavfile - The file path of an IDL save file produced by
;                mage_fitstd which contains the fit to the standard
;   objstr     - The object structure generated from the mage_script
;                extraction routines
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
;   REJ - Set the rejection threshold for fit points.  If the value of
;         an order fit at a particular wavelength is less than this
;         percentage of the maximum of the fit to the order then that
;         point will be flagged.  
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mage_flux,'gd108cal.sav',obj_strct,rej=0.02
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
;   x_calcfit
;
; REVISION HISTORY:
;   16-Jun-2008 CLW

pro mage_flux,calsavfile,objstr,rejfrac=rejfrac,CHK=CHK
  
  if not keyword_set(rejfrac) then rejfrac=0.05
  restore,file=calsavfile
  norders=n_elements(objstr)
  if n_elements(tot_fit) ne norders then $
     message,"Wrong number of orders in calfile"
  
  for ord=0,norders-1 do begin
     gdwv = where(objstr[ord].wave > 1000.0)
     ofit = x_calcfit(objstr[ord].wave[gdwv],fitstr=(tot_fit[ord]))
     nfit = n_elements(ofit)
     ;;;;;;;;;
     ;; It has occurred that ofit is greater than zero in the wings
     ;; (as the Legendre polynomial shoots up) which messes with the
     ;; fluxing and order coaddition (and all steps afterwards).
     ;; So going to only consider the ofit where continuous and not
     ;; curving upwards (hopefully this leaves enough of an overlap)
     ;; KLC  18 Jan 2011
     ;;;;;;;;;
     gdfit = where(ofit gt 0.0)
     wvc = 61750./(6.+(norders-1.-ord)) ; center of order (see mage_1dspec)
     if 6+(norders-1.-ord) le 8 then begin
        ;; This approximation breaks down towards the red orders (ord ~
        ;; 14); so going to check on the real data
        gd = where(objstr[ord].fx[gdwv] gt 0.)
        if gd[0] eq -1 then begin
           print,'mage_flux: no good data in this order ',6.+(norders-1.-ord)
           gd = gdfit           ; use it all
        endif 
        mx = max(ofit[gdwv[gd]],imx)
        wvc2 = objstr[ord].wave[gdwv[gd[imx]]]
        if abs(wvc2-wvc) gt 100. then wvc = wvc2 ; Trust data
        ;; Kludge for Order 14 :: JXP 16Aug2012
        if ord EQ 14 then wvc = 9800.  ;; Not good if you really have flux beyond 1micron! 
     endif 
     ;if ord EQ 1 then stop
     dwv = wvc/(6.+(norders-1.-ord))*1.5             ; sigma outwards
     tmp = min(objstr[ord].wave[gdwv[gdfit]]-wvc,icent,/absolute) ; peak
     rng = gdfit[icent] + lindgen(nfit-gdfit[icent]) ; right-hand side
     bd = where(ofit[rng] le 0. or objstr[ord].wave[gdwv[rng]] gt wvc+dwv)
     if bd[0] ne -1 then $
        tmp = min(ofit[rng[0:bd[0]]],ihi) $ ; find the dip
     else tmp = min(ofit[rng],ihi)
     hi = rng[ihi]               ; last acceptable
     ;if ord EQ 1 then stop
     rng = lindgen(gdfit[icent]) ; left-hand side
     bd = where(ofit[rng] le 0. or objstr[ord].wave[gdwv[rng]] lt wvc-dwv,nbd) 
     if bd[0] ne -1 then begin
        tmp = min(ofit[rng[bd[nbd-1]:gdfit[icent]-1]],ilo)
        ilo = bd[nbd-1]+ilo  
     endif else tmp = min(ofit[rng],ilo)
     lo = rng[ilo]              ; first acceptable
     ;if ord EQ 1 then stop

     ;; Visual Check
     ;if ord EQ 14 then $
      ;x_splot,objstr[ord].wave[gdwv],ofit,xtwo=objstr[ord].wave[gdwv[lo:hi]],ytwo=ofit[lo:hi],/block
;     if 6.+(norders-1.-ord) gt 7 and 6.+(norders-1.-ord) lt 9 then stop
     ;stop
     ofit[0:lo] = 0.0           ; lo, hi might be negative; exclude
     ofit[hi:nfit-1] = 0.0
     ;;;;;;;;;
      izero=WHERE(ofit GT 0.0,nzero)
     imid=izero[nzero/2]
     indx=lindgen(nfit)
     xcen=long_find_nminima(-ofit,indx,nfind=5,width=10L,minsep=5L $
                            ,ypeak=ypeak,npeak=npeak)
     icen_good=where(xcen GT 0.1*izero[0] AND xcen LT 0.9*izero[nzero-1L],ncen)
     IF ncen EQ 0 THEN message,'Problem with your sensitivity function' $
     ELSE BEGIN
        xcen=xcen[icen_good] 
        ypeak=ypeak[icen_good]
        omax=max(-ypeak,kmax)
        imax=xcen[kmax]
     ENDELSE
     ileft=0
     FOR kk=imax,0,-1L DO BEGIN
        IF ofit[kk] LT rejfrac*omax THEN BEGIN
           ileft=kk
           BREAK
        ENDIF
     ENDFOR
     iright=nfit-1L
     FOR kk=imax,nfit-1L,1L DO BEGIN
        IF ofit[kk] LT rejfrac*omax THEN BEGIN
           iright=kk
           BREAK
        ENDIF
     ENDFOR
     qgood=objstr[ord].wave gt 1000 AND ofit GT 0.0 AND $
           indx GE ileft AND indx LE iright
     ;; FIX THIS LATER WITH A DERIVATIVE TO FIND THE MAX
     ;;qgood=objstr[ord].wave gt 1000 AND ofit GT 0.0 
     igd=WHERE(qgood,ngd)
     objstr[ord].flux[igd]=objstr[ord].fx[igd]/objstr[ord].exp/ofit[igd]
     objstr[ord].sky[igd]=objstr[ord].sky[igd]/objstr[ord].exp/ofit[igd]
     objstr[ord].sig[igd]=sqrt(objstr[ord].var[igd])/objstr[ord].exp/ofit[igd]
     objstr[ord].nosig[igd]= $
        sqrt(objstr[ord].novar[igd])/objstr[ord].exp/ofit[igd]
     ;; mask bad points
     ibd=WHERE(qgood EQ 0 OR objstr[ord].var LE 0.0,nbd)
     IF nbd GT 0 THEN BEGIN
        objstr[ord].flux[ibd]=0.0
        objstr[ord].sky[ibd]=0.0
        objstr[ord].sig[ibd]=0.0
        objstr[ord].nosig[ibd]=0.0
     ENDIF
     IF KEYWORD_SET(CHK) THEN BEGIN
        iwv=where(objstr[ord].WAVE GT 1000.0)
        min_wv=min(objstr[ord].wave[iwv])
        max_wv=max(objstr[ord].wave[iwv])
        clr=getcolor(/load)
        plot,objstr[ord].wave,ofit,yr=[0.0,omax],xrange=[min_wv,max_wv] 
        oplot,objstr[ord].wave[igd],ofit[igd],color=clr.red
        wait,3
     ENDIF
  ENDFOR
END
