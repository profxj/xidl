;; Return the indices of specified frames from a mage structure
FUNCTION MAGE_FRAMES, mage, frames

  igd = WHERE(frames GE 0, nframes)
  ikeep = lonarr(nframes)
  FOR iframe = 0L, nframes-1L DO $
     ikeep[iframe] = WHERE(mage.FRAME EQ frames[igd[iframe]])
  RETURN, ikeep
END
