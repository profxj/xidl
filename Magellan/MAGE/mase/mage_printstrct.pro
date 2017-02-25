pro mage_printstrct, mage

  science = mage[where(mage.exptype EQ "SCIENCE" and mage.obj_id GE 0)]

  objids = uniq(science.obj_id)

  for id=0, n_elements(objids)-1 do begin

     print, science[objids[id]].obj_id, " ", science[objids[id]].object

  endfor

end
