;; Converts ion numbers into a label
function x_lblion, ion, flag, LATEX=latex

  ;;
  if not keyword_set(flag) then flag = 0

  cup = '!u'
  cdown = '!N'
  if keyword_set(LATEX) then begin
     cup = '^{'
     cdown = '}'
  endif
  
  ;;
  getabnd, elm, ion[0], abnd, flag=1
  case ion[1] of
     1: addc = cup+'0'+cdown
     2: addc = cup+'+'+cdown
     3: addc = cup+'++'+cdown
     4: addc = cup+'3+'+cdown
     5: addc = cup+'4+'+cdown
     else: stop
  endcase
  case ion[1] of
     1: ionc = 'I'
     2: ionc = 'II'
     3: ionc = 'III'
     4: ionc = 'IV'
     else: stop
  endcase

  lbl = elm+addc
  lbl2 = elm+ionc

  case flag of
     0: ulbl = lbl
     1: ulbl = lbl2
     else: stop
  endcase

  return, ulbl
end
