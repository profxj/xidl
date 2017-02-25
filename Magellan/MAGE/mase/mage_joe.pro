files=['Objstr0080.fits','Objstr0082.fits','Objstr0257.fits','Objstr0259.fits']

  FOR j=0L,n_elements(files)-1L DO BEGIN
     obj_str= mrdfits(files[j],1) 
     IF j EQ 0 THEN allframes=obj_str $
     ELSE allframes = [allframes,obj_str]
  ENDFOR
  
END

mage_combspec, allframes, fspec

END
