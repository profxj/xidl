function xlc, lin
  nlin = strlen(lin)
  i = strpos(lin,' ')
  if(i EQ -1) then return, nlin else return, i
end
