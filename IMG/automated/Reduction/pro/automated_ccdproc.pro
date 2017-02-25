;+
; NAME:
;   automated_ccdproc
;
; PURPOSE:
;   Generate fully reduced images (bias,dark,gain,flat)
;
;    
; CALLING SEQUENCE:
; 
; automated_ccdproc, plan, SKYFIT=, GZIP=, CLOBBER=
;
; INPUTS:
;   
; plan - the structure with relevant information
; 
; 
; OPTIONAL INPUTS:
;  
; skyfit       - subtract a second order polynomial (not active)
; gzip         - zip the final data
; clobber      - overwrite data that have been already reduced
;
; OUTPUTS: 
;
; Creates a bunch of reduced science frames. It also creates an
; additional structure that keeps track of what has been done.  
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; 
;
; EXAMPLES:
; 
; Rereduce data with the structure STR
; automated_ccdproc, STR, /CLOBBER, 
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;    05-Jun-2012  Written and revised by MF
;
;-
;------------------------------------------------------------------------------
;; 


PRO automated_ccdproc, plan, skyfit=skyfit, gzip=gzip, clobber=clobber


  ;;open the plan structure 
  str=mrdfits(plan,1,/silent)
    
  ;;find science, bias and flat
  bias=where(str.type EQ 'bia', nbias)
  dark=where(str.type EQ 'dar', ndark)
  flat=where(str.type EQ 'fla', nflats)
  science=where(str.type EQ 'ima', nsci)
  

  ;;if clobber
  if keyword_set(clobber) then spawn, 'rm -f '+plan+'status.str'

  ;;check if status str exists
  aa=FILE_SEARCH(plan+"status.str",count=cc)
  IF (cc EQ 0 ) THEN BEGIN
      ;;create status structure 
      splog, "Creating new STATUS structure"
      status={BIASMED:0,DARK:[0.,0.],FLAT:0,FILFLAT:["","","","","",""],IMG:MAKE_ARRAY(N_ELEMENTS(str.name),/integer,value=0)}
      mwrfits, status, plan+"status.str", /create
  ENDIF ELSE BEGIN
      splog, "STATUS structure found. Continuing from there..."
      status=mrdfits(plan+"status.str",1,/silent)
  ENDELSE
  



;;----------------
;;make bias
;;----------------
 
  ;;check if bias already exists 
  
  IF(status.biasmed EQ 1) THEN  BEGIN
     splog, "Using exisiting median bias."
     biasmedian=mrdfits(plan+"medbias.fits",0,/silent)
  ENDIF   
  
  IF(status.biasmed EQ 2) THEN  splog, "No bias frame. Use overscan only!"


  IF(nbias GT 0 AND status.biasmed EQ 0) THEN BEGIN
     automated_makebias, str, biasmedian, bias, PLAN=plan
     splog, "Done with median bias"
     status.BIASMED=1
  ENDIF 
  
  IF(nbias EQ 0 AND status.biasmed EQ 0) THEN BEGIN
     splog, "No bias frame. Use overscan only"
     status.BIASMED=2
  ENDIF



  ;;update structure
  mwrfits, status, plan+"status.str", /create

  ;;----------
  ;;make dark
  ;;---------

  

  if(status.dark[0] ne 0.) then  splog, "Using exisiting dark."
  if(status.dark[0] eq -9999) then  splog, "No darks found...I can live with that for now!"


  if(ndark gt 0 and status.dark[0] eq 0) then begin
      if(status.biasmed eq 1) then automated_makedark, str, dark, status=status, bias=biasmedian
      if(status.biasmed eq 2) then automated_makedark, str, dark, status=status
      splog, "All done with darks"
  endif 
  
  if(ndark eq 0 and status.dark[0] eq 0) then begin
     splog, "No darks found..."
     status.dark[0]=-9999
  endif
 
 

  ;;update structure
  mwrfits, status, plan+"status.str", /create

  
  ;;---------------
  ;;make flat 
  ;;---------------


  ;;set min and max values to accept a flat as good
  minmax=[1000.,41000.]
  

  ;;check if flats already exists 
  if(status.flat eq 1) then  splog, "Using exisiting flats."
  if(status.flat eq 2) then begin  ;;this parts needs to be tested more (MF June 2012)
     splog, "No flats found... Try with archive!"
     automated_makeflats, str, flat, status=status, plan=plan, minmax=minmax
     status.flat=1
 endif

  
  if(nflats gt 0 and status.flat eq 0) then begin
      if(status.biasmed eq 1) then automated_makeflats, str, flat, status=status, bias=biasmedian, plan=plan,minmax=minmax
      if(status.biasmed eq 2) then automated_makeflats, str, flat, status=status, plan=plan, minmax=minmax
      splog, "All done with flats"
      status.flat=1
  endif 
  
  if(nflats eq 0 and status.flat eq 0) then begin
     splog, "No flats found... I cannot reduce data without flats!"
     splog, "Exiting.. Nothing will happen!"
     status.flat=2
     return
  endif
  
  ;;update structure
  mwrfits, status, plan+"status.str", /create
  
  
  ;;--------------------
  ;;make science frame
  ;;--------------------
  
  if(nsci gt 0) then begin
      
      ;;loop on each science frame
      for img=0, nsci-1 do begin
      
          if(status.img[img] eq 1) then  splog, "Found previously reduced frame ", str.name[science[img]] else begin
              
              if(status.biasmed eq 1) then automated_makescience, str, science[img],status=status,plan=plan, $
                      bias=biasmedian, skyfit=skyfit, gzip=gzip
              if(status.biasmed eq 2) then automated_makescience, str, science[img],status=status,plan=plan, $
                      skyfit=skyfit, gzip=gzip
              
              ;;update structure
              status.img[img]=1
              mwrfits, status, plan+"status.str", /create
          endelse
          
      endfor
      
      
  endif else begin
      print, "No science frame!"
      return
  endelse
  
  
  splog, "All done with data reduction! "
  
  
end
