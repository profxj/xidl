;+ 
; NAME:
; wfccd_combfspec
;    Version 1.0
;
; PURPOSE:
;   Converts SDSS eigenfunctions into WFCCD
;
; CALLING SEQUENCE:
;   
;   wfccd_combfspec, wffspec
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   XSIZE      - Size of gui in screen x-pixels (default = 1000)
;   YSIZE      - Size of gui in screen y-pixels (default = 600)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_combfspec, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   16-May-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_combfspec, fspec_fil, out_fil, PARSE=parse, OBJDIR=objdir

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_combfspec, fspec_fil, out_fil [v1.0]'
    return
  endif 

; Optional Keywords

; Open fspec files
  nfspec = n_elements([fspec_fil])
  if nfspec EQ 1 then begin
      print, 'wfccd_combfspec: Assuming this is a wildcard!'
      a = findfile(fspec_fil, count=cnt)
      if cnt EQ 1 then begin
          print, 'wfccd_combfspec: Assuming this is a wildcard!'
          return
      endif else begin
          fspec_fil=a
          nfspec = cnt
      endelse
  end

  for q=0L, nfspec-1 do begin
      wfccd_wrfspec, wffspec, fspec_fil[q], /read
      if q EQ 0 then tot_fspec = wffspec else tot_fspec=[tot_fspec,wffspec]
  endfor

; Sort on obj name
  list = strarr(n_elements(tot_fspec))
  for q=0L,n_elements(tot_fspec)-1 do $
    list[q] = strtrim(tot_fspec[q].slit_id,2)+strtrim(tot_fspec[q].obj_id,2)

; Unique test
  
  uobj = uniq(list, sort(list))
  nuobj = n_elements(uobj)
  if nuobj NE n_elements(list) then begin
      print, 'wfccd_combfspec: Not all obj are unique!'
      stop
  endif

; Parse z
  if keyword_set( PARSE ) then begin
      sz_prs = size(parse, /dimensions)
      for i=0L,sz_prs[0]-1 do begin
          obj = strtrim(parse[i,0],2)+x_objnumid(parse[i,1])
          indx = where(list EQ obj, nindx)
          if nindx EQ 0 then stop
          case parse[i,2] of
              -1: tot_fspec[indx].zans.z = -1.
              else: stop
          endcase
      endfor
  endif

; Write File
  ;; Sort
  srt = sort(tot_fspec.slit_id)
  wfccd_wrfspec, tot_fspec[srt], out_fil

; OBJ Dir
  if keyword_set( OBJDIR ) then begin
      for i=0L,n_elements(tot_fspec)-1 do begin
          a = where(strlen(strtrim(tot_fspec[i].obj_fil,2)) NE 0,na)
          if na NE 0 then tot_fspec[i].obj_fil[a] = $
            OBJDIR+tot_fspec[i].obj_fil[a]
      endfor
  endif
  
  
; Write File
  wfccd_wrfspec, tot_fspec[srt], out_fil

  print, 'wfccd_combfspec: All done!'
  return
end


