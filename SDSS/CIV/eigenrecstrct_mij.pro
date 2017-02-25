;+
; NAME:
;   EIGENRECSTRCT_MIJ
;     Version 2.0
;
; AUTHOR:
;   Melodie M. Kao
;   Massachusetts Institute of Technology
;   Kavli Institute for Astrophysics and Space Research
;   77 Massachusetts Avenue, Building 37-287, Cambridge, MA 02139 
;   melodie.kao@alum.mit.edu
;   mkao@caltech.edu
;
; PURPOSE:
;   Fits predetermined eigenspectra to a spectrum. Primary function:
;   finds m_ij matrix
;
;
; CALLING SEQUENCE:
;   m_invert = EIGENRECSTRCT_MIJ( weight, allEigflux, [m_ij=, fail=,
;                                 /silent]) 
;
;
; DESCRIPTION: 
;
;   EIGENRECONSTRUCT will shift spectra into a given rest frame and
;   fit given eigenspectra to the spectra.  The eigenspectra must be
;   an ascii table with the first column as z=0.0 wavelengths and the
;   second column as fluxes (example: See Yip C. W. et al., 2004, AJ,
;   128, 2603 for an example of correctly formatted files).  They will
;   be fed into the program as a list of file names that include the
;   complete path to the eigenspectra.
;   The spectra will be returned as 10^(lambda_0 + i*dLoglambda_0).  
;
;
; INPUTS:
;
;   weight    --  weighting used in fit
;   allEigflux -- Array of eigenspectra to be used in fitting
;
;
; RETURNS:
;
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;   /silent -- don't print anything
;
;   fail      -- Integer for each spectrum fitted:
;                -1 = Small pivot element used to invert matrix,
;                     significant accuracy probably lost
;                 0 = Matrix inversion successful, accuracy pretty good
;                 1 = Singular array, inversion is invalid, no
;                 continuum
;   m_ij= -- some intermediate matrix
;
; COMMENTS: 
;   Written to be used with eigenreconstruct.pro.
;
; EXAMPLES:
;   See usage in eigenreconstruct.pro
;
; PROCEDURES/FUNCTIONS CALLED:
;   interpspec.pro              FIRE/Spextool/pro/
;
; MODIFICATION HISTORY:
;  25 April 2011   Written by M. Kao
;
;-
; Copyright (C) 2011, Melodie Kao
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-
;------------------------------------------------------------------------------

function eigenrecstrct_mij_rdfil, eigspecfil, outfile=outfile
  ;; Just read in everything
  ;; eigspecfil can either be a string array of file names or a single
  ;; name that will be read into a string array
  if n_params() ne 1 then begin
     print,'Syntax - eigenrecstrct_mij_rdfil( eigspecfil, [outfile=])'
     return,-1
  endif 

  neig = (SIZE(eigspecfil,/DIM))[0]
  if neig eq 1 then begin
     ;; Assume have to read list
     readcol,eigspecfil,eigname,format='a',/silent
     neig = (size(eigname,/dim))[0]
  endif else eigname = eigspecfil

  neigbins = FILE_LINES(eigname[0]) ; assume all files the same (should be 2582)
  eigenArr = DBLARR(neigbins,neig+1,/nozero) ; could be dangeraus if didn't knwo files same
  
  if keyword_set(outfile) then begin
     kywd = 'EIGEN'+strtrim(lindgen(neig),2)
     prs = strsplit(eigname[0],'/',count=nprs,/extract)
  endif 

  FOR i = 0L, neig-1 DO BEGIN
     OPENR, unitj, eigname[i], /GET_LUN ; this faster than readcol
     FOR j = 0, neigbins-1 DO BEGIN
        READF, unitj, wv, fx
        if i eq 0 then eigenArr[j,0] = wv
        eigenArr[j,i+1] = fx
     ENDFOR                     ; loop j=neigbins
     FREE_LUN, unitj

     if keyword_set(outfile) then begin
        prs = strsplit(eigname[i],'/',count=nprs,/extract)
        sxaddpar,header,kywd[i],prs[nprs-1],'Included eigenspectrum' 
     endif 
  ENDFOR                        ; loop i=neig

  if keyword_set(outfile) then begin
     mwrfits,eigenArr,outfile,header,/create,/silent
     print,'eigenrecstrct_mij_rdfil(): created ',outfile
  endif 

  return, eigenArr

end                             ; eigrenrecstrct_mij_rdfil()


function eigenrecstrct_mij_interpspec, eigenArrFull, wave, modes=modes, $
                                       zQSO=zQSO, zin=zin, list=list, fits=fits
 ;; Optional Keywords
 ;;  zin       -- Redshift of frame that spectrum should be shifted to.
 ;;              Keep zframe = 0.0 if you want the spectrum to be
 ;;              shifted to rest frame of the spectrum object. Function
 ;;              works by first shifting to restframe of the quasar and
 ;;              then shifting to actual requested restframe. 
 ;;              Default = 0.0
 ;;  zQSO= -- redshift used to shift the input wavelength (wave)
 ;;  modes     -- Integer array of eigenmodes that you want to
 ;;               use. Eigenmode order is determined by order of text
 ;;               list fed into function (eigspecfil); if list order is:
 ;;               [eigfile1, eigfile4, eigfile2, eigfile0, eigfile3] and
 ;;               you give and array of [0,2,3], then the modes used
 ;;               will be [eigfile1, eigfile2, eigfile0].  Default is
 ;;               all of the supplied modes.
 ;;

;  if n_params() ne 2 then begin
;     print,'Syntax - eigenrecstrct_mij_interpspec(eigenArrFull, wave, [modes=, zQSO=, zin=])'
;     return,-1
;  endif 
  ;; I ask for two different z values because first eigenrecstrct_mij
  ;; will shift to the observed frame of the qso and then to the
  ;; restframe of the LRG.  If the spectra are already in their restframe
  ;; or are only going to be stacked on top of each other, then first
  ;; shift the spectra to the restframe, then feed into
  ;; eigenrecstrct.pro zin = zQSO = 0.0. 

  ;; Truncate
  if keyword_set(modes) then begin
     if (size(modes,/dim))[0] eq 0 then $ ; Not array so use [0:modes+1]
        eigenArr = eigenArrFull[*,0:modes] $  ; instead of modes-1
     else eigenArr = eigenArrFull[*,[0,1+modes]] ; assume array
  endif else eigenArr = eigenArrFull

  npix = (size(wave,/dim))[0]
  neig = (size(eigenArr,/dim))[1]   ; [1st dimen = wave, 2nd dimen = eigspec + 1]
  allEigflux  = DBLARR(neig-1,npix,/nozero)

  IF keyword_set(zin) THEN BEGIN      
     if not keyword_set(zQSO) then $
        stop,'eigenrecstrct_mij_interpsec(): zin and zQSO must be set together.'
     ;; shift to observed frame of qso then shift to restframe of object
     eigenArr[*,0] = eigenArr[*,0]*(1.0+zQSO)/(1.0+zin)
  ENDIF   
  
  FOR i = 1L, neig-1 DO BEGIN     ; b/c [*,0] is wavelength
     INTERPSPEC, eigenArr[*,0], eigenArr[*,i], wave, eigflux
     allEigflux[i-1,*] = eigflux
  ENDFOR                        ; loop i=neig
  return, allEigflux

end                             ; eigenrecstrct_mij_interspec()


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function eigenrecstrct_mij, weight, allEigflux, m_ij=m_ij, $
                            FAIL=fail, SILENT=silent
;  if n_params() ne 2 then begin
;     print,'Syntax - eigenrecstrct_mij(weight, allEigflux, [m_ij=, FAIL=, /SILENT])'
;     return,-1
;  endif 

  sz = size(allEigflux,/dim)
  neig = sz[0]
  npix = sz[1]
  
  ;; Set arrays
  warr = replicate(1.d,neig) # weight ; blow up to be [neig,npix]
  m_ij = matrix_multiply(alleigflux,alleigflux*warr,/btranspose) ; m_ij = [neig, neig]

  m_invert = INVERT(m_ij,status,/DOUBLE)
  case status of
     0: BEGIN
        if not keyword_set(silent) then $
           print, "eigenrecstrct_mij: successful matrix inversion"
        fail = 0
     END
     1: BEGIN
        if not keyword_set(silent) then $
           print, "eigenrecstrct_mij: singular array; inversion is invalid"
        fail = 1
     END
     2:  BEGIN
        if not keyword_set(silent) then $
           print, "eigenrecstrct_mij: small pivot element used; significant accuracy probably lost"
        fail = -1
     END
     else: stop,'eigenrecstrct_mij: unknown status ',status
  endcase                       ; status
  
  return, m_invert
end
