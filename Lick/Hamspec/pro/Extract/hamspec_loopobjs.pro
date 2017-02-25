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
;   /USEARCHFLAT -- Use the archival dewar 4 pixel flat, currently
;   getenv('XIDL_DIR') + '/Lick/Hamspec/calibs/hamspec_dewar4_MF_2013may_1x1.fits.gz'
;
; COMMENTS:
;
; EXAMPLES:
;  hamspec_loopobjs,strct,setup
;
; PROCEDURES/FUNCTIONS CALLED:
;  hamspec_proc
;  hamspec_fntobj
;  hamspec_extract
;
; REVISION HISTORY:
;   01-May-2013 Written by BPH
;-
;------------------------------------------------------------------------------
pro hamspec_loopobjs,hamspec_strct,setup,FLATFIL=flatfil,CLOBBER=clobber,USEARCHFLAT=usearchflat

  ;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hamspec_loopobjs, hamspec_strct,setup'
      return
  endif 

  if not keyword_set(CLOBBER) then clobber=0
  if not keyword_set(USEARCHFLAT) then usearchflat=0
  if not keyword_set(FLATFIL) then flatfil=''

   if keyword_set(USEARCHFLAT) and not keyword_set(FLATFIL) then flatfil = getenv('XIDL_DIR') + '/Lick/Hamspec/calibs/hamspec_dewar4_MF_2013may_1x1.fits.gz'

  obj_ind = where(hamspec_strct.setup eq setup and hamspec_strct.obj_id gt 0 ,nobj)
  
  if nobj eq 0 then begin
     print, "No object files with a setup of "+strtrim(setup,2)
     return
  endif

  inds = uniq(hamspec_strct[obj_ind].obj_id,sort(hamspec_strct[obj_ind].obj_id))
  ninds = n_elements(inds)
  if ninds eq 1 then inds = hamspec_strct[setup_ind[0]].obj_id

  
  for cur_objind=0,ninds-1 do begin
     cur_obj = hamspec_strct[obj_ind[inds[cur_objind]]].obj_id
     hamspec_proc,hamspec_strct,setup=setup,obj=cur_obj,FLATFIL=flatfil,CLOBBER=clobber
     hamspec_fntobj, hamspec_strct, setup, cur_obj, objaper = [1., 1.], CLOBBER=clobber
     hamspec_extract, hamspec_strct, setup, cur_obj,  /reschk, /bbox, CLOBBER=clobber
  endfor

end
