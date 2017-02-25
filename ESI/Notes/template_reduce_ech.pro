pro reduce_ech_all, esi, slit, OSTRT=ostrt, CUAR=cuar, NOOSCAN=nooscan, ORDRS=ordrs

 if abs(slit-0.75) LT 1e-3 then FLUXFIL = 'Extract/sens_esi0079.idl' else stop
 if not keyword_set(OSTRT) then OSTRT = 0L
 if not keyword_set(FILSTD) then filstd = 'Extract/Obj_esi0013.fits'

 ;; Find all unique object numbers
 gd = where(esi.type EQ 'OBJ' and abs(esi.slit-slit) LT 0.01 $
            and esi.flg_anly NE 0 and esi.mode EQ 2, ngd)
 if ngd EQ 0 then return
 uni_objid = esi[gd[uniq(esi[gd].obj_id, sort(esi[gd].obj_id))]].obj_id
 nobjid = n_elements(uni_objid)

 ;; Loop
 for qq=OSTRT,nobjid-1 do begin

     obj_id = uni_objid[qq]
     print, 'reduce_all: Working on obj_id = ', obj_id
     idx = where(esi[gd].obj_id EQ obj_id and abs(esi[gd].slit-slit) LT 0.01 AND $
                esi[gd].mode EQ 2 and esi[gd].flg_anly NE 0)
     objind = gd[idx]
     print, 'reduce_all: Working on obj = ', esi[objind[0]].obj

     ;; Proc image
     esi_echproc, esi, objind, /CLOBBER, NOOSCAN=nooscan

     ;; Zap CRs
     if n_elements(objind) GT 1 then esi_echobjcr, esi, obj_id

     ;; Find object
     esi_echfndobj, esi, obj_id, chk = chk, nfind = 1, /USESTD, FWHM = 5.0, FILSTD=filstd

     ;; Sky subtract
     esi_echskysub, esi, obj_id, fchk = chk, bordr = 3;, CUAR=cuar

     ;; Refine trace
     esi_echfndobj, esi, obj_id, chk = chk, nfind = 1, $
                    /USESTD, /SKYSUB, FWHM=5.0, FILSTD=filstd

     ;; Trace and extract object.
     esi_echextobj, esi, obj_id, chk=fchk, /OPTIMAL, ORDRS=ordrs

     ;; Now combine, flux, and coadd spectra from this night
     if keyword_set(ORDRS) then eordr = [ordrs[0], ordrs[n_elements(ordrs)-1]]
     esi_echcombspec, esi, obj_id, obj_nm = 'a', ordrs=eordr
     esi_echfluxfin, esi, obj_id, fluxfil = fluxfil, obj_nm = 'a'
     esi_echcoaddfin, esi, obj_id, obj_nm = 'a', /SKY, /NOVAR, ORDRS=eordr
 endfor

 return
end

