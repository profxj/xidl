;+ 
; NAME:
; hamspec_loopobjs   
;     Version 1.1
;
; PURPOSE:
;
; This simply loops through all of the objects observed with a given
; setup. It processes each object and does the one-d extraction.
; 
; CALLING SEQUENCE:
;   
;  hamspec_loopobjs,struct, setup
;  
;
; INPUTS:
;   [hamspecstruct] - hamspec structure array
;   [setup]   -  Setup identifier 
;
; RETURNS:
;  Structure, setup 
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /CLOBBER -- stomp on existing reductions
;   FLATFIL -- Flat field (pixel flat) for image processing
;
; COMMENTS:
;
; EXAMPLES:


;
; PROCEDURES/FUNCTIONS CALLED:
;  hamspec_strct
;  hamspec_setup
;  hamspec_ar()
;  hamspec_wrstrct
;
;  hamspec_loopobjs
;
; REVISION HISTORY:
;   01-Jun-2013 Written by BPH
;-
;------------------------------------------------------------------------------


pro hamspec_pipeline, indir, CLOBBER=clobber, USEBIAS=usebias, USESETUP=usesetup, STRETCH=stretch

  ;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'hamspec_pipeline, indir'
      return
  endif 



  if not keyword_set(USESETUP) then usesetup=0
  if not keyword_set(CLOBBER) then clobber=0
  if not keyword_set(USEBIAS) then usebias=0
  if not keyword_set(STRETCH) then stretch=0

  if usesetup then begin
     hamspec = hamspec_ar()
     hamspec_setup,hamspec
  endif else begin
     ; need to make one
     hamspec_strct,hamspec
     hamspec_setup,hamspec
     hamspec_wrstrct,hamspec
  endelse
; need to get list of setups


  if hamspec[0].ccd ne 'e2vCCD203-824kx4kthin' then begin
     print,'Pipeline currently only works with Dewar 4'
     return
  endif

  minsetup = min(hamspec.setup)
  maxsetup = max(hamspec.setup LT 50)

  if maxsetup LT minsetup then begin
     print, 'No data setups, only calibration'
     return
  endif

; make the pixel flat
; find the files identified as MFLT 
  calibration_ind = where(hamspec.flg_anly NE 0 AND $
                          strtrim(hamspec.type,2) EQ 'MFLT',nmflt )

  if nmflt gt 0 then begin
; find all unique setups with MFLT

     mflt_setups=uniq(hamspec[calibration_ind].setup,sort(hamspec[calibration_ind].setup))
     nmflt_setups = n_elements(mflt_setups)
; only one setup, likely 50
     if nmflt_setups eq 1 then begin 
        first = calibration_ind[0]
        mflt_setup = hamspec[first].setup
        hamspec_mkmflat,hamspec,mflt_setup, CLOBBER=clobber, USEBIAS=usebias, $
                        SMOOTH=smooth
     endif else begin
        ; this really shouldn't happen but
        mflt_setup = mflt_setups[0]
        for ii=0,nmflt_setups-1 do begin
           hamspec_mkmflat,hamspec,mflt_setups[ii], CLOBBER=clobber, USEBIAS=usebias, $
                           SMOOTH=smooth
        endfor
     endelse
; read in the output of hamspec_mkmflat, this will be passed to
; hamspec_proc at a later step

     flatfil_name = 'Flats/Flat_'+strtrim(mflt_setup,2)+'_M.fits' 
  endif else begin
; no MFLT files,  use archived file
     flatfil_name = getenv('XIDL_DIR') + '/Lick/Hamspec/calibs/hamspec_dewar4_MF_2013may_1x1.fits.gz'


  endelse 

; all setups
  all_setup_inds = uniq(hamspec.setup,sort(hamspec.setup))
; remove setups that are calibrations, >50
  setup_inds = where(hamspec[all_setup_inds].setup lt 50,nsetup)
  for ii=0,nsetup-1 do begin
     setup = hamspec[all_setup_inds[setup_inds[ii]]].setup
     outfil = hamspec_getfil('qtz_fil', setup, /name, CHKFIL=chkf)     
     if chkf eq 0 then begin
        tflt_ind = where(hamspec.setup EQ setup $
                         and hamspec.type EQ 'TFLT', $
                         ntflt)
        if ntflt gt 0 then begin 
           hamspec_mkflats,hamspec,setup,'TFLT', USEBIAS=usebias, $
                           AVG=avg, CLOBBER=clobber, _EXTRA=extra
           hamspec_edgeflat, hamspec, setup, CLOBBER=clobber
           hamspec_nrmflat, hamspec, setup, CLOBBER=clobber
        endif 
        
     endif 
     hamspec_mkarc, hamspec, setup
     ll =  getenv('XIDL_DIR')+'/Spec/Arcs/Lists/ham_dewar4_thar.lst' 
     hamspec_fitarc,hamspec,setup,linlist=ll,STRETCH=stretch, CLOBBER=clobber
     hamspec_fit2darc, hamspec, setup, nocoeff=3, CLOBBER=clobber
     hamspec_tracearc, hamspec, setup, CLOBBER=clobber
     hamspec_fittrcarc, hamspec, setup, CLOBBER=clobber
     hamspec_mkaimg, hamspec, setup,  CLOBBER=clobber
     hamspec_slitflat, hamspec, setup,  CLOBBER=clobber

     hamspec_loopobjs, hamspec, setup, FLATFIL=flatfil_name ; loop through the objects in the setup

  endfor


end
  


