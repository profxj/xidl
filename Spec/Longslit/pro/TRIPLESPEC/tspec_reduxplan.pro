;;  Reduces a series of images from a plan file
;;   Can do Telluric's too
PRO tspec_reduxplan, plan_file, telluric=telluric, CLOBBER=clobber

  ;; Read in plan file
  readcol, plan_file, dum1, filenm, flg, type, name, exp, dum2, $
           dum3, dum4, dum5, dum6, airmass, idx1, idx2, $
           format='A,A,I,A,A,F,F,F,I,A,I,F,I,I' 
  nframe = n_elements(type)

   ;; Telluric vs. Object
   if keyword_set(TELLURIC) then otype = 'telluric' else otype = 'object'

   ;; Strip quotes off of type as need be
   for ii=0L,nframe-1 do begin
      ipos = strpos(type[ii],'"')
      if ipos LT 0 then continue
      type[ii] = strsplit(type[ii], '"', /extr)
   endfor

   ;; Find all of the objects or tellurics
   if Not keyword_set(telluric) then $
      obj = where(strmatch(type,otype)) $
   Else  obj = where(strmatch(type,otype))

    ;; Find unique object names
   uni_nm = name[ obj[uniq( name[obj], sort(name[obj]))]] 
   nuni = n_elements(uni_nm)

  ;; Loop
  for qq=0L,nuni-1 do begin
    ;; Find all with that name
    gd = where( strmatch(name[obj],uni_nm[qq]) )
    print, 'Working on '+otype,  uni_nm[qq]
    ;; 
    mxi1 = max( idx1[obj[gd]] )
    for jj=0L,mxi1 do begin
       ;;
       gd2 = where( strmatch(name[obj],uni_nm[qq]) and idx1[obj] EQ jj)
       ;; 
       mxi2 = max( idx2[obj[gd2]] )
       case mxi2 of
          0: begin
             if jj GT 0 then begin
                ;; Use the previous frame
                gd3 = where( strmatch(name[obj],uni_nm[qq]) and idx1[obj] EQ (jj-1))
                gd2 = [gd2,max(gd3)]
             endif else begin
                print, 'tspec_reduxplan: Only 1 file.  This will not work!!'
                stop
             endelse
          end
          3: 
          else: print, 'tspec_reduxplan: Was expecting 4 files but found only', mxi2+1
       endcase
       for ii=0L,mxi2,2 do begin
          ;; Pair of filenames
          if ii+1 GT mxi2 then begin
             ;; Odd number.  Will use previous frame and hope for ABBA
             if ii EQ 0 then begin
                if n_elements(gd2) EQ 2 then rdx_file = filenm[ obj[gd2[ii+lindgen(2)]]] $
                else stop
             endif else rdx_file = filenm[ obj[gd2[ii-lindgen(2)]]]
             ;stop ;; Not well tested -- JXP on October 1, 2013
          endif else rdx_file = filenm[ obj[gd2[ii+lindgen(2)]]]
          ;; Reduce
          print, 'Reducing files ', rdx_file
          ;; Telluric?
          if keyword_set(TELLURIC) then begin
             itel = obj[gd2[ii]]
             robj = where(strmatch(type,'object'))
             mn = min(abs(itel-robj), imn)
             oidx = robj[imn]
             if oidx LT itel then opair = [oidx-1,oidx] else opair = [oidx,oidx+1]
             obj_fil  = filenm[opair]
             ;; Reduce first
             tspec_redux_pair, obj_fil, telluric=rdx_file, clobber=clobber
             ;; Reduce second
             tspec_redux_pair, obj_fil, telluric=reverse(rdx_file), clobber=clobber
          endif else begin
             ;; First one
             print, 'tspec_reduxplan: Reducing first exposure'
             tspec_redux_pair, rdx_file
             ;; Second one
             print, 'tspec_reduxplan: Reducing second exposure'
             tspec_redux_pair, reverse(rdx_file)
          endelse
          ;stop
       endfor
    endfor
  endfor

  return
end
