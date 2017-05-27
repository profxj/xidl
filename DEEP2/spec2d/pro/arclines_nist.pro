; 2002-Jul-23 - DPF
; read spectral line lists from NIST database
;
;All lines obtained from
;
;NIST Atomic Spectra Database
;http://physics.nist.gov
;
;http://physics.nist.gov/cgi-bin/AtData/main_asd

;Indexes of refraction are derived for expressions given by 
;  E.R. Peck and K. Reeder, J. Opt. Soc. Am. 62, 958 (1972).
; 
; Note that I am NOT converting back to vacuum properly here. 
;
function nist_read, fname

  str5 = strarr(5)
  spawn, 'wc -l '+fname, result
  nline = (strsplit(result[0], ' ', /extract))[0]-6
  str = strarr(nline)
  openr, rlun, fname, /get
  readf, rlun, str5
  readf, rlun, str
  free_lun, rlun
  w = where(strmid(str, 0, 2) ne '  ', nline)
  str = str[w]
  col = strsplit(str[0], '|')
  element = strmid(str, col[0], col[1]-1)
  lambdastr = strmid(str, col[1], col[2]-col[1]-1)
  lambda = double(lambdastr)
  intensstr = strmid(str, col[2], col[3]-col[2]-1)
  intensity = float(stregex(intensstr,'[0-9]+',/extract))
  intnote = (stregex(intensstr,'[0-9]+([^0-9])',/extract,/sub))[1, *]
  intnote = transpose(intnote)

  line1 = {lambda: 0.0d, intensity:0.0, element:'', ionization:'', $
           note:'', quality:''}

  line = replicate(line1, nline)
  line.element = strmid(element, 0, 2)
  line.intensity = intensity
  line.lambda = lambda
  line.ionization = strcompress(strmid(element, 2), /remove)
  line.note = intnote

  return, line
end

pro arclines_nist, a, mindist=mindist, minint=minint

  mindist = 4
  minint = 5

; -------- set path
  deep_dir = getenv('DEEP_DIR')+'/'
  if deep_dir eq '/' then message, 'you must set $DEEP_DIR'
  path = deep_dir+'spec2d/etc/nist_data/'

; Do not include He anymore
;  a = nist_read(path+'He.dat')
  arcfiles=['Ne.dat','Ar.dat','Kr.dat','Xe.dat']
  arcfiles=['Cd.dat','Hg.dat','Ne.dat','Ar.dat','Kr.dat','Xe.dat']


  narc=n_elements(arcfiles)

  for i=0,narc-1 do begin
      atmp=nist_read(path+arcfiles[i])
      bright = where(atmp.intensity GE minint, na)
      atmp = atmp[bright]

      ind = sort(atmp.lambda)
      atmp = atmp[ind]
      atmp.quality = 'GOOD'
      lam = atmp.lambda
      lam0 = shift(lam, -1)
      lam1 = shift(lam, +1)
      for j=0, na-1 do begin
         min = min(abs([lam0[j], lam1[j]]-lam[j]))

        if min LT mindist then atmp[j].quality = 'BLEND'
        if strcompress(atmp[j].note, /remove) ne '' then atmp[j].quality = 'BAD'
     endfor 

      if i eq 0 then a=atmp else a=[a,atmp]
 endfor
;  a = nist_read(path+'Ne.dat')
;  a = [a, nist_read(path+'Ar.dat')]
;  a = [a, nist_read(path+'Kr.dat')]
;  a = [a, nist_read(path+'Xe.dat')]

;  bright = where(a.intensity GE minint, na)
;  a = a[bright];

;  ind = sort(a.lambda)
;  a = a[ind]
;  a.quality = 'GOOD'
;  lam = a.lambda
;  lam0 = shift(a.lambda, -1)
;  lam1 = shift(a.lambda, +1)
;  for i=0, na-1 do begin
;     min = min(abs([lam0[i], lam1[i]]-lam[i]));

;     if min LT mindist then a[i].quality = 'BLEND'
;     if strcompress(a[i].note, /remove) ne '' then a[i].quality = 'BAD'
;  endfor 
;    n = 1.000275d

  ind = sort(a.lambda)
  a = a[ind]

  airlambda=a.lambda
  vactoair,airlambda 
  na=n_elements(airlambda)

  for i=0, na-1 do print, airlambda[i], a[i].intensity, a[i].quality, $
    a[i].element, format='(F10.4,I6,A6,A3)'

  listexists=findfile(deep_dir+'spec2d/etc/lamp_NIST_blue.dat') NE ''
  
  if listexists eq 0 then begin
      openw,2,deep_dir+'spec2d/etc/lamp_NIST_blue.dat'
      for i=0, na-1 do printf, 2,airlambda[i], a[i].intensity, a[i].quality, $
        a[i].element, format='(F10.4,I6,A6,A3)'
      close,2
  endif
  return
end


