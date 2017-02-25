;+
; NAME:
;   mmt_fixhead
;
; PURPOSE:
;  Fix a minor 'bug' in the MMT BCS header
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:

; REVISION HISTORY:
;   20-Apr-2005  Written by J. Hennawi Berkeley
;-
;------------------------------------------------------------------------------

pro mmt_fixhead, root, grating, instrum

  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mmt_fixhead, root, grating, [instrum], [v1.0]'
      return
  endif 

  if not keyword_set(INSTRUM) then instrum = 'mmtbluechan'

  ;; Grating check
  grtng = ['800GPM', '300_5300']
  mt = where(grating EQ grtng,nmt)
  if nmt EQ 0 then begin
      print, 'Not a valid grating! ', grating
      return
  end

  fil = findfile(root+'*', count=nfil)

  for qq=0L,nfil-1 do begin
      ;; Read
      dat = xmrdfits(fil[qq],0,head,/sile)

      ;; Update
      sxaddpar, head, 'DISPERSE', grating
      sxaddpar, head, 'INSTRUME', instrum
      
      ;; Write
       slen = strlen(fil[qq])
       flg_gz = 0
       if strmid(fil[qq],slen-3) EQ '.gz' then $
         outfil = strmid(fil[qq],0,slen-3) else outfil = fil[qq]
      mwrfits, dat, outfil, head, /create
      spawn, 'gzip -f '+outfil

      print, 'mmt_fixhead: Updated ', outfil
  endfor
  
  return
END
