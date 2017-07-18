;+
; NAME:
;   SDSS_IONCHK
;     Version 1.0
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
;   Checks to see if there are other potential ion absorption line
;   matches in a system if given a redshift of a detected ion.
;
;
; CALLING SEQUENCE:
;   sdss_ionchk, civfil,dblt_name=dblt_name,NEWCIVFIL=newcivfil,
;                DVTOL=dvtol, DEBUG=debug
;
; DESCRIPTION:
;   Ions used for matching (in order of importance--if no more space
;   remains in the structure to store all matches, the matches are
;   stored in order of importance): 
;               1548.195        CIV
;               1550.770        CIV
;               2796.352        MgII
;               2803.531        MgII
;               1393.755        SiIV
;               1402.770        SiIV
;               1334.5323       CII
;               1526.7066       SiII
;               2600.1729       FeII
;               2586.6500       FeII
;               2382.765        FeII 
;               2344.214        FeII 
;               2374.4612       FeII 
;               1260.4221       SiII
;               1302.1685       OI
;               1608.4511       FeII
;               1670.7874       AlII
;               1854.7164       AlIII
;               1862.7895       AlIII
;
; INPUTS:
;   civfil -- FIT file of all CIV candidate systems.  Structure should
;             be defined using sdsscivstrct__define.pro
;
;
; RETURNS: 
;
;
; OUTPUTS:
;   A FIT file (replacing original input file) of all input systems, 
;   with additional information stored in obj.wrest and obj.zabs_orig. 
;   The .wrest tag stores the rest wavelengths of ions that matched 
;   detected centroids, USING INDICES [2:10].  The .zabs_orig tag
;   stores the z of the detected ion-centroid matches, USING INDICES [2:10].
;
;
; OPTIONAL KEYWORDS:
;   ionlist -- Enumerate ions of interest to measure EW
;              (limit). (Default: long ion list.)
;
;   dblt_name -- Default ion used as reference point to other
;                ions being searched.  Default CIV
;
;   newcivfil -- File name of FIT file storing original input
;                information and up to 8 detected ion-centroid
;                matches.  If given, original input file will not be
;                overwritten.
;
;   dvtol     -- Velocity tolerance.  Any centroids found within EITHER
;                the ion rest wavelength +/- a wavelength interval
;                defined by this velocity tolerance OR the equivalent
;                width wave limits (therefore nonsymmetric tolerances
;                are ok) are accepted as a match. (Default: 150 km/s)
;
;   /DEBUG    -- For each ion-centroid match found, prints out:
;                  Index #      
;                  Ion #
;                  ion_rest    [angstroms]   
;                  z_ion 
;                  z_dblt
;                  ion_offset  [km/s]
;                  blue_tol    [km/s]
;                  red_tol     [km/s]
;
; OPTIONAL OUTPUTS:
;    A new FIT file (name defined by optional input newcivfil)
;    of all input systems, with additional information
;    stored in obj.wrest and obj.zabs_orig.  The .wrest tag stores the
;    rest wavelengths of ions that matched detected centroids, USING
;    INDICES [2:10].  The .zabs_orig tag stores the z of the detected
;    ion-centroid matches, USING INDICES [2:10].
;
; COMMENTS: 
;
;
; EXAMPLES:
;  
;
; PROCEDURES/FUNCTIONS CALLED:
;  
;
; MODIFICATION HISTORY:
;  15 Jul 2011   Written by M. Kao
;  13 Sep 2011   Revamped a bit, KLC
;   7 Jul 2014   Updated, KLC
;
;-
; Copyright (C) 2011, Melodie Kao
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-
;------------------------------------------------------------------------------



pro sdss_ionchk, civfil,dblt_name=dblt_name,NEWCIVFIL=newcivfil,    $
                 DVTOL=dvtol, use_cflg=use_cflg, DEBUG=debug, ionlist=ionlist
  
  if  N_params() LT 1  then begin 
     print,'Syntax - sdss_ionchk, civfil, [dblt_name=,newcivfil=,dvtol=,'
     print,'                      use_cflg=, /debug]'
     return
  endif 
  sdssdir = sdss_getsdssdir()

  ;; Set initial values
  cinv = 1./2.998e5
  IF NOT KEYWORD_SET(dvtol) THEN dvtol = 150.0 ; km/s
  IF NOT KEYWORD_SET(newcivfil) THEN BEGIN
     print, 'sdss_ionchk: Overwrite existing file?: y(es), n(o).'
     done = 0
     while not done do begin
        ch = get_kbrd(1)        ; w/(1): wait for character
        case byte(strlowcase(ch)) of
           110: begin           ; n
              done = 1
              print,  'Please specify an output file as NEWCIVFIL= '
              return
           end
           121: done = 1        ; y
        endcase
     ENDWHILE
     outcivfil = civfil
  ENDIF else outcivfil = newcivfil
  IF NOT KEYWORD_SET(dblt_name) then dblt_name = 'CIV'
  if size(dblt_name,/type) eq 8 then dblt = dblt_name $
  else dblt = dblt_retrieve(dblt_name)

  
  if not keyword_set(ionlist) then begin
     ;; Define ion order of importance/interest
     ionlist = [1548.195 , $    ; CIV
                1550.770 , $    ; CIV
                2796.352 , $    ; MgII
                2803.531 , $    ; MgII
                1393.755 , $    ; SiIV
                1402.770 , $    ; SiIV
                1334.5323, $    ; CII
                1526.7066, $    ; SiII
                2600.1729, $    ; FeII
                2586.6500, $    ; FeII
                2382.765 , $    ; FeII
                2344.214 , $    ; FeII 
                2374.4612, $    ; FeII 
                1260.4221, $    ; SiII
                1302.1685, $    ; OI
                1608.4511, $    ; FeII
                1670.7874, $    ; AlII
                1854.7164, $    ; AlIII
                1862.7895, $    ; AlIII
                2852.9633]      ; MgI
     
     ;; Exclude doublet of interest
     gd = where(ionlist ne dblt.wvI and ionlist ne dblt.wvII,ngd)
     ionlist = ionlist[gd]
  endif
  ionsize = (SIZE(ionlist, /DIM))[0]
  
  if size(civfil,/type) eq 7 then $
     civstr = xmrdfits(civfil,1, /SILENT) $
  else civstr = civfil
  civstr.abslin_fil = strtrim(civstr.abslin_fil,2)
  ncivstr = (SIZE(civstr, /DIM))[0]
  nionmax = (size(civstr.wrest,/dim))[0]

  ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ;; %%%%%%%%%%%%%%%%%%%%%%%%%% DEBUG INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ;; %%%%%%                                                            %%%%%%
  IF KEYWORD_SET(debug) THEN BEGIN 
     print, ''
     print, 'i'+dblt.ion,'iIon','ion_rest','z_ion','z_'+dblt.ion,'dvabs',  $
            'dv_blue','dv_red', $
            format='(a5,1x,a4,2x,a9,2(1x,a7),1x,3(a7,1x))'
  ENDIF
  ;; %%%%%%                                                            %%%%%%
  ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  prev_abslin_fil = ''          ; don't read in more often than necessary
  FOR icivstr = 0L, ncivstr-1 DO BEGIN

     ;; Check whether need to read in new abslin structure
     if prev_abslin_fil ne civstr[icivstr].abslin_fil then begin
        test      = FILE_SEARCH( sdssdir+civstr[icivstr].abslin_fil+'*', $
                                 COUNT=ntest )
        IF ntest EQ 0 THEN STOP, "sdss_ionchk: --FILE NOT FOUND--", $
                                 sdssdir+civstr[icivstr].abslin_fil

        ;; Get all absorption centroids found by sdss_fndlin.pro
        abslin = xmrdfits( test[0], 1, /SILENT) 
        if keyword_set(use_cflg) then begin
           if abslin.cflg ne use_cflg then begin
              print,'sdss_ionchk: EW and WVLIM not saved since cflg does not match'
              cflg_mismtch = 1
           endif else cflg_mismtch = 0
           abslin.cflg = use_cflg
        endif else cflg_mismtch = 0
        cindx = fix(alog(abslin.cflg)/alog(2))
        ;; Save for next interation check
        prev_abslin_fil = civstr[icivstr].abslin_fil
     endif 
     

     ;; Figure out which wrest index to start at
     iflg = WHERE( civstr[icivstr].wrest EQ 0.00, nflg)
     IF nflg EQ 0 THEN begin
        print, "sdss_ionchk: No space to store match information: ",$
               civstr[icivstr].qso_name,$
               string("zabs=",civstr[icivstr].zabs_orig[0],format='(a,f7.5)')
        continue                ; EXIT
     endif 


     iflg = iflg[0]
     IF civstr[icivstr].zabs_final[0] GT 0. THEN                        $
        zabs = civstr[icivstr].zabs_final[0]                              $
     ELSE zabs = civstr[icivstr].zabs_orig[0]
     ctrdrest = abslin.centroid[0:abslin.ncent[cindx]-1,cindx] / (1.0 + zabs)
     
     iion = 0
     while iion lt ionsize and iflg lt nionmax do begin
        ;; Test if ion already found and saved
        match = where(abs(civstr[icivstr].wrest-ionlist[iion]) lt 1.e-4)
        if match[0] ne -1 then begin
           if keyword_set(debug) then $
              print,'sdss_ionchk debug: ion already in structure ',ionlist[iion]
           iion++
           continue
        endif 

        ;; Define wavelength tolerances. Use a lower bound and a
        ;; higher bound (blue_tol and red_tol, respectively).
        ;; Define each as the greater of either the wavelength
        ;; limits that determine the EW or dvtol
        IF civstr[icivstr].wvlim_final[0,0] GT 0. THEN BEGIN
           blue_tol = ABS(dblt.wvI*(1.0+zabs)-civstr[icivstr].wvlim_final[0,0])
           red_tol  = ABS(dblt.wvI*(1.0+zabs)-civstr[icivstr].wvlim_final[0,1])
        ENDIF ELSE BEGIN    $
           blue_tol = ABS(dblt.wvI*(1.0+zabs)-civstr[icivstr].wvlim_orig[0,0])
           red_tol  = ABS(dblt.wvI*(1.0+zabs)-civstr[icivstr].wvlim_orig[0,1])
        END
        default_tol = ionlist[iion]*cinv*dvtol
        blue_tol = default_tol > blue_tol
        red_tol  = default_tol > red_tol
        
        ;; Find where potential ion matches are detected
        match = WHERE( ctrdrest-ionlist[iion] GE -1.0*blue_tol AND        $
                       ctrdrest-ionlist[iion] LE red_tol     , nmatch)
        IF nmatch NE 0 THEN BEGIN
           ;; Even if there are multiple matches, we only need one
           ;; match.  Select the abslin centroid closest to the
           ;; ion wavelength.
           civstr[icivstr].wrest[iflg]     =  ionlist[iion]
           bestmatch = MIN( ABS(ctrdrest[match]-ionlist[iion]), imatch)
           match = match[imatch]
           civstr[icivstr].zabs_orig[iflg] = $
              abslin.centroid[match,cindx] / ionlist[iion] - 1

           if cflg_mismtch eq 0 then begin
              ;; Able to save this information because input structure
              ;; and desired use_cflg match
              ;; So store and turn in to rest EW
              cnst = 1./(1. + civstr[icivstr].zabs_orig[iflg])
              civstr[icivstr].ew_orig[iflg] = abslin.ew_orig[match]*cnst
              civstr[icivstr].sigew_orig[iflg] = abslin.sigew_orig[match]*cnst
              civstr[icivstr].wvlim_orig[iflg,*] = abslin.wvlim_orig[match,*]
              civstr[icivstr].ewflg[iflg] = sdss_getewflg(/custom) ; 64 to indicate questionable
           endif 
           
           
           ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           ;; %%%%%%%%%%%%%%%%%%%%% DEBUG INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           ;; %%%%%%                                                  %%%%%%
           IF KEYWORD_SET(debug) THEN BEGIN
              print, icivstr, iflg, civstr[icivstr].wrest[iflg], $
                     civstr[icivstr].zabs_orig[iflg], zabs, $
                     (civstr[icivstr].zabs_orig[iflg] - zabs)/(1+zabs)*2.998e5,$
;                     (ctrdrest[match]-ionlist[iion])*2.998e5 / $
;                     abslin.centroid[match,cindx], $
                     -1.0*blue_tol*2.998e5 / abslin.centroid[match,cindx], $
                     red_tol* 2.998e5 / abslin.centroid[match,cindx], $
                     format='(i5,1x,i4,2x,f9.4,2(1x,f7.5),1x,3(f7.2,1x))'
              
              if iflg+1 eq nionmax and iion+1 ne ionsize then $
                 print,'sdss_ionchk debug: out of room but not done with ions'
              
           ENDIF 
           ;; %%%%%%                                                  %%%%%%
           ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           ;; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
           
           ;; Increment
           iflg++
        ENDIF                   ; If ion matches are found

        iion++
     ENDWHILE                   ; Loop through ions

  ENDFOR                        ; Looping through all of the candidates
  
  
  if size(outcivfil,/type) eq 7 then begin
     mwrfits, civstr, outcivfil, /create, /silent
     spawn,'gzip -f '+outcivfil
     print,'sdss_ionchk: wrote to ',outcivfil
  endif else begin
     if keyword_set(newcivfil) then newcivfil = civstr $
     else civfil = civstr
  endelse
  
end
