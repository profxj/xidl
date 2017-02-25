function m_mkobjstr, norders

  norders=15
  tmp = {objstrct}
  objstr = replicate(tmp, norders)
  objstr.UT = ' '
  objstr.field = ' '
  objstr.img_fil = ' '
  objstr.arc_fil = ' '
  objstr.exp = 0
  objstr.field = ' '
  objstr.flg_anly = 1

  return, objstr

END
