;+ 
; NAME:
; esi_reduce_all_obj   
;     Version 1.0
;
; PURPOSE:
;    Reduce all esi frames after all processing of Bias, Arcs, Flats,
;    and STD star is finished. 
;
; CALLING SEQUENCE:
;   
;  esi_reduce_all_obj, esi, slit, [obj_id], [OSTRT], [/chk], [/fchk], [nfind]
;
; INPUTS:
;   esi      -  ESI structure
;   slit     -  The slit size, (e.g. 0.75)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   [FLUXFIL=]  - Specify the File created to flux the spectrum.
;   [obj_id=]   - Object ID  (e.g. 0L, 1L, etc) (Does NOT take an array)
;   [OSTRT=]    - Start from a later object, based on this number of
;                 objects. (probably same as obj_id). Keeps going after.
;   /FCHK       - Turns on the check for Trace and Extract object -
;                 Recommend doing this at least the first time
;   /CHK        - Show many checks on many steps (do if there is a problem)
;   [nfind=]    - Specify number of objects to find in the
;                 slit. (e.g. 1L, 2L, 3L) 
;                 Recommend doing this only with obj_id set!
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
; esi_reduce_all_obj, esi, slit, fluxfil=fluxfil, /fchk
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   ??-???-2009 reduce_all_obj Written by JXP
;   08-Jan-2010 renamed and modified by MR
;
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_reduce_all_obj, esi, slit, obj_id=obj_id, OSTRT=ostrt,FLUXFIL=fluxfil, $
                        chk=chk, fchk=fchk, nfind=nfind

if keyword_set(nfind) and not keyword_set(obj_id) then begin
   print, 'Are you sure you want to find multiple objects in a slit for all object IDs?'
   print, 'Press .c to continue'
   stop
endif

if keyword_set(fchk) then fchk=1 else fchk=0
if keyword_set(chk) then chk=1 else chk=0
if not keyword_set(nfind) then nfind=1L
if not keyword_set(FLUXFIL) then FLUXFIL = 'Extract/sens_esi0039.idl'
if not keyword_set(OSTRT) then OSTRT = 0L


 ;; Find all unique object numbers
 gd = where(esi.type EQ 'OBJ' and abs(esi.slit-slit) LT 0.01 $
            and esi.flg_anly NE 0 and esi.mode EQ 2, ngd)
 if ngd EQ 0 then return
 uni_objid = esi[gd[uniq(esi[gd].obj_id, sort(esi[gd].obj_id))]].obj_id
 nobjid = n_elements(uni_objid)

; Makes the next loop only go over the obj_id specified.
 if keyword_set(obj_id) then begin
    OSTRT = obj_id
    nobjid = obj_id+1
 endif
 
 ;; Loop
 for qq=OSTRT,nobjid-1 do begin
    
    obj_id = uni_objid[qq]
    print, 'reduce_all: Working on obj_id = ', obj_id
    idx = where(esi[gd].obj_id EQ obj_id and abs(esi.slit-slit) LT 0.01 AND $
                esi[gd].mode EQ 2 and esi[gd].flg_anly NE 0)
    objind = gd[idx]
    print, 'reduce_all: Working on obj = ', esi[objind[0]].obj
    
      ;; Proc image
     esi_echproc, esi, objind, /CLOBBER

     ;; Zap CRs
     if n_elements(objind) GT 1 then esi_echobjcr, esi, obj_id

     ;; Find object
     esi_echfndobj, esi, obj_id, chk = chk, nfind = nfind, /USESTD, FWHM = 5.0

     ;; Sky subtract
     esi_echskysub, esi, obj_id, fchk = chk, bordr = 3

     ;; Refine trace
     esi_echfndobj, esi, obj_id, chk = chk, nfind = nfind, /USESTD, /SKYSUB, FWHM = 5.0

     ;; Trace and extract object.
     esi_echextobj, esi, obj_id, chk=fchk, /OPTIMAL

     ;; Now comibine, flux, and coadd spectra from this night
     esi_echcombspec, esi, obj_id, obj_nm = 'a'
     esi_echfluxfin, esi, obj_id, fluxfil = fluxfil, obj_nm = 'a'
     esi_echcoaddfin, esi, obj_id, obj_nm = 'a', /SKY, /NOVAR
     
     if nfind eq 2 then begin
        esi_echcombspec, esi, obj_id, obj_nm = 'b'
        esi_echfluxfin, esi, obj_id, fluxfil = fluxfil, obj_nm = 'b'
        esi_echcoaddfin, esi, obj_id, obj_nm = 'b', /SKY, /NOVAR
     endif 
     
     if nfind eq 3 then begin
        esi_echcombspec, esi, obj_id, obj_nm = 'c'
        esi_echfluxfin, esi, obj_id, fluxfil = fluxfil, obj_nm = 'c'
        esi_echcoaddfin, esi, obj_id, obj_nm = 'c', /SKY, /NOVAR
     endif 
     
     endfor
 
 return
end

