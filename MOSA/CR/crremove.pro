;+ 
; NAME:
; crremove
;    Version 1.1
;
; PURPOSE:
;  Loops on the various CCDs in the MOSA image and use la_cosmic
;  to identify the CRs
;
; CALLING SEQUENCE:
;  crremove, list, msklist, iraffil, NITER=niter, SONLY=sonly
;
; INPUTS:
;  list -- List of input images
;  msklist -- List of mask names to modify
;  iraffil -- Not sure this is necessary
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  NITER=  -- Number of iterations for la_cosmic
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Apr-2003 Written by JXP
;------------------------------------------------------------------------------
pro crremove, list, msklist, iraffil, NITER=niter, SONLY=sonly

  ;; Read list
  readcol, list, filnm, FORMAT='A'
  readcol, msklist, mskfil, FORMAT='A'

  if not keyword_set( NITER ) then niter = 3L

  ;; IRAFFIL
  if keyword_set(iraffil) then begin
      close, /all
      openw, 1, iraffil
  endif
  ;; Loop on list
  nfil = n_elements(filnm)

  for qq=0L,nfil-1 do begin
      ;; Loop on sub-images
      for i=1L,8 do begin
          ;; Read image
;          img = xmrdfits(strtrim(filnm[qq],2)+'.fits', i, /silent)

          ;; Run la_cosmic
          msknm = strtrim(mskfil[qq],2)+'_'+strtrim(i,2)
          if not keyword_set( SONLY ) then begin
              la_cosmic, [strtrim(filnm[qq],2)+'.fits'], $
                masklist=[msknm+'.fits'], $
                gain=3.1, $
                readn=5.6, sigclip=4.5, sigfrac=0.5, objlim=2., niter=niter, $
                blocksize=1024L, indx=i
          endif
          ;; IRAFFIL (convert to pl)
          printf, 1, 'imcopy '+msknm+'.fits '+msknm+'.pl'
          printf, 1, 'imdel '+msknm+'.fits'
      endfor
      ;; More IRAF
      printf, 1, 'crplusbpmask "', strmid(filnm[qq],3), '"'
      for i=1L,8 do $
        printf, 1, '#fixpix ', strtrim(filnm[qq],2)+'['+strtrim(i,2)+'] BPM'
  endfor

  ;;
  close, /all

  return
end
