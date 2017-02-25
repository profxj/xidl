FUNCTION LUCI_FILENO, filename, prefix = prefix
  
  split1 = strsplit(filename, '\.fits*', /extract, /regex)
  temp = split1[0]
  temp1 = strsplit(temp, '\.', /extract, /regex)
  indx = temp1[n_elements(temp1)-1]
  ;split1 = strsplit(filename, 'fits*', /extract)
  ;ifits = WHERE(strmatch(split1, 'fits*'))
  ;indx = split1[ifits[0]-1L]

  split2 = strsplit(filename, '/', /extract)
  nsplit = n_elements(split2)
  prefix = strsplit(split2[nsplit-1L], '\.fits', /extract, /regex)

  RETURN, indx
END
