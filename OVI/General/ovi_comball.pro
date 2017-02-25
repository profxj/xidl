;+ 
; NAME:
; ovi_comball
;    Version 1.1
;
; PURPOSE:
;   Takes a list of files representing all of the observations for
; a given field and combines them to make one galaxy summary file.
; Unless /NOUNIQ is set, the routine will check that multiple observations
; of the same object gave the same redshift and complain if not the case.
;
; CALLING SEQUENCE:
; ovi_comball, fspec_fil, obj_fil, out_fil, PARSE=, NOUNIQ=, ZSET=
;
; INPUTS:
;  fspec_fil  -- List of WFCCD spectra files (string array)
;  obj_fil    -- Object files (string array)
;
; RETURNS:
;
; OUTPUTS:
;  out_fil  -- Combine WFCCD spectra file
;
; OPTIONAL KEYWORDS:
;  PARSE   -- 2-element long array which manually sets the flags of 
;             specific objects (mainly to reject bad redshifts)
;  /NOUNIQ -- Skip the unique checks.  
;  ZSET=   -- List of [name, redshift] to manually set the redshift 
;             of an obj
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   ovi_comball, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro ovi_comball, fspec_fil, obj_fil, out_fil, PARSE=parse, $
                 NOUNIQ=nouniq, ZSET=zset


;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
      'wfccd_combfspec, fspec_fil, obj_fil, out_fil, PARSE=, ZSET= [v1.0]'
    return
  endif 

; Directories
  a = findfile('GSpec/..', count=count)
  if count EQ 0 then file_mkdir, 'GSpec'
  a = findfile('GSpec/Slits/..', count=count)
  if count EQ 0 then file_mkdir, 'GSpec/Slits'
  a = findfile('GSpec/Final/..', count=count)
  if count EQ 0 then file_mkdir, 'GSpec/Final'
  a = findfile('GSpec/Objfil/..', count=count)
  if count EQ 0 then file_mkdir, 'GSpec/Objfil'

; Optional Keywords
  if not keyword_set( INSTR ) then instr = 'W'

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

  ;; LOOP
  for q=0L, nfspec-1 do begin
      ;; Open fspec
      wfccd_wrfspec, wffspec, fspec_fil[q], /read
      ;; Truncated objfil
      trunc_obj = obj_fil[q,*,0]
      all_ipos = strpos(obj_fil[q,*,0], 'Extract')
      for i=0L,n_elements(all_ipos)-1 do begin
          if all_ipos[i] NE -1 then $
            trunc_obj[i] = strmid(obj_fil[q,i,0],all_ipos[i])
      endfor
      ;; Loop on Spec
      for jj=0L,n_elements(wffspec)-1 do begin
          a = where(strlen(strtrim(wffspec[jj].obj_fil,2)) NE 0,na)
          for ii=0L,na-1 do begin
              ;; Grab the date
              indx = where(wffspec[jj].obj_fil[a[ii]] EQ trunc_obj,nind)
              if nind EQ 0 then stop else indx = indx[0]
              ;; Parse out the frame number and make new file names
              ipos = strpos(obj_fil[q,indx,0],'Obj')
              newobjfil = 'GSpec/Objfil/O'+strtrim(obj_fil[q,indx,1],2)+ $
                '_'+instr+strmid(obj_fil[q,indx,0],ipos+4,3)+'.fits'
              newfinfil = 'GSpec/Final/F'+strtrim(obj_fil[q,indx,1],2)+ $
                '_'+instr+strmid(obj_fil[q,indx,0],ipos+4,3)+'.fits'
              oldfinfil = strmid(obj_fil[q,indx,0],0,all_ipos[indx])+ $
                'Final/f_ccd'+strmid(obj_fil[q,indx,0],ipos+4,3)+'.fits.gz'
              if keyword_set( SVOBJFIL ) then begin
                  b = where(newobjfil EQ svobjfil,nb)
              endif else nb = 0
              ;; Update Obj and Write!
              if nb EQ 0 then begin
                  if not keyword_set( SVOBJFIL ) then svobjfil = newobjfil $
                  else svobjfil = [svobjfil, newobjfil]
                  ;; Read obj structure
                  objstrct = xmrdfits(obj_fil[q,indx,0],1,$
                                      structyp='specobjstrct',/silent)
                  ;; Reset spec2d_fil
                  objstrct.spec2d_fil = newfinfil
                  ;; Reset Slit file
                  oldslit = objstrct[0].slit_fil
                  ipos = strpos(oldslit, 'ts_')
                  newslit = 'GSpec/Slits/S'+strtrim(obj_fil[q,indx,1],2)+ $
                    '_'+instr+strmid(oldslit,ipos+3,2)+'.fits'
                  objstrct.slit_fil = newslit
                  ;; Write New Obj file
                  print, 'ovi_comball: Creating ', newobjfil
                  mwrfits, objstrct, newobjfil, /create
                  spawn, 'gzip -f '+newobjfil
                  ;; Copy over Final file
                  print, 'ovi_comball: Copying over ', oldfinfil
                  spawn, 'cp '+oldfinfil+' '+newfinfil+'.gz'
              endif
              ;; Reset name
              wffspec[jj].obj_fil[a[ii]] = newobjfil
          endfor
      endfor
      ;; Update Date
      if q EQ 0 then tot_fspec = wffspec else tot_fspec=[tot_fspec,wffspec]
  endfor

; Sort on obj name
  list = strarr(n_elements(tot_fspec))
  for q=0L,n_elements(tot_fspec)-1 do $
    list[q] = strtrim(tot_fspec[q].slit_id,2)+strtrim(tot_fspec[q].obj_id,2)

; Unique test
  
  uobj = uniq(list, sort(list))
  if not keyword_set( NOUNIQ ) then begin
      nuobj = n_elements(uobj)
      if nuobj NE n_elements(list) then begin
          print, 'wfccd_combfspec: Not all obj are unique!'
          stop
      endif
  endif else begin
      msk = intarr(n_elements(list))
      msk[uobj] = 1B
      ;; Loop on all obj
      for i=0L,n_elements(list)-1 do begin
          ;; Skip unique
          if msk[i] EQ 1 or msk[i] EQ -1 then continue
          ;; Check parse
          if keyword_set( PARSE ) then begin
              a = where(tot_fspec[i].slit_id EQ parse[*,0], na)
              if na NE 0 then stop  ; It is probably safe to continue
          endif
          ;; Check z values
          a = where(tot_fspec.slit_id EQ tot_fspec[i].slit_id AND $
                    tot_fspec.obj_id EQ tot_fspec[i].obj_id, na)
          if na LE 1 then stop
          for j=1L,na-1 do begin
              diff = abs(tot_fspec[a[j]].zans.z - tot_fspec[i].zans.z)
              if diff GT 0.001 then stop
          endfor
          ;; Print
          print, 'ovi_comball: Multiple measurements for', list[i]
          ;; Take first element
          msk[i] = 1
          msk[a[1:*]] = -1
      endfor
      gd = where(msk EQ 1)
      list = list[gd]
      tot_fspec = tot_fspec[gd]
  endelse

  ;; Parse z
  if keyword_set( PARSE ) then begin
      sz_prs = size(parse, /dimensions)
      for i=0L,sz_prs[0]-1 do begin
          obj = strtrim(parse[i,0],2)+x_objnumid(parse[i,1])
          indx = where(list EQ obj, nindx)
          if nindx EQ 0 then stop
          case parse[i,2] of
              -1: tot_fspec[indx].zans.z = -1.  ;; Confused blend
              -2: begin  ;; Bad slit
                  tot_fspec[indx].zans.z = -1.
                  tot_fspec[indx].flg_anly = 1
              end
              -3: begin  ;; Star
                  tot_fspec[indx].zans.z = 0.
              end
              else: stop
          endcase
      endfor
  endif

  ;; Zset
  if keyword_set( ZSET ) then begin
      for i=0L, n_elements(zset)-1 do begin
          indx = where(list EQ zset[i].name, nindx)
          if nindx EQ 0 then stop else $
            tot_fspec[indx].zans.z = zset[i].z
      endfor
  endif

; Write File
  ;; Sort
  srt = sort(tot_fspec.slit_id)
  wfccd_wrfspec, tot_fspec[srt], out_fil

; Write File
  wfccd_wrfspec, tot_fspec[srt], out_fil

  print, 'ovi_comball: All done!'
  return
end


