function mage_findobj_old, sciimg, skyimg, scivar, orderfil=orderfil, hdr=hdr, QA=qafil, std=std, nopca=nopca, chk=chk

  IF NOT FILE_TEST('QA', /DIR) THEN FILE_MKDIR, 'QA'

   if (not keyword_set(qafil)) then begin
     
      qafil = "QA/QAFindObj.ps"
   ENDIF

   ordr_strct = mrdfits(orderfil, 1)
   obj_strct = m_mkobjstr(15)
   obj_strct.exp = sxpar(hdr,"EXPTIME")

   if (keyword_set(STD) OR keyword_set(NOPCA))then begin
      m_fntobj, obj_strct, ordr_strct, sciimg-skyimg, scivar, qafil, std=std
   endif else begin
      m_fntobj, obj_strct, ordr_strct, sciimg-skyimg, scivar, qafil
   endelse

   ;xatv, (sciimg-skyimg) * (skyimg GT 0), min=-20, max=100, /block
   ;for qq=0L,14 do xatvplot, obj_strct[qq].trace[0:2048-1], findgen(2048)
   
   ;Suggested by Rob, added by JJB
   IF KEYWORD_SET(CHK) THEN BEGIN
   tmp = sciimg-skyimg
   for qq=0L, 14 do begin
      tmp[obj_strct[qq].trace[0:2048-1], findgen(2048)] = -1000
   endfor
   xatv, tmp, /block, min=-100, max=100
ENDIF




   return, obj_strct

end
