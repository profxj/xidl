FUNCTION jfh_fit_reject, x, y, ncoeff, inmask = inmask, outmask = outmask $
                         , ivar = ivar1, maxrej = maxrej $
                         , TOL = TOL, YFIT = YFIT, SIGREJ = SIGREJ $
                         , FUNC = FUNC, XMIN = XMIN, XMAX = XMAX

IF NOT KEYWORD_SET(SIGREJ) THEN SIGREJ = 3.0
IF NOT KEYWORD_SET(TOL) THEN TOL = 1.0e10
nx = n_elements(x)
IF KEYWORD_SET(INMASK) THEN OUTMASK = INMASK $
ELSE outmask = lonarr(nx)+1
use_inds = where(outmask, nuse)
IF NOT KEYWORD_SET(MAXREJ) THEN maxrej = nuse/2
; Now fit a polynomial to the fwhm and reject outliers
IF nuse GT 1 THEN BEGIN
    FOR j = 0, maxrej-1L DO BEGIN
        IF KEYWORD_SET(ivar1) THEN ivar = ivar1[use_inds]
        xy2traceset, x[use_inds], y[use_inds], xset, invvar = ivar $
                     , ncoeff = ncoeff, maxiter = 0, /silent, FUNC = FUNC $
                     , XMIN = XMIN, XMAX = XMAX $
                     , yfit = yfit
        resids = y[use_inds]-yfit
        sig = djsig(resids)
        rej_inds = WHERE(abs(resids) GE sigrej*sig OR abs(resids) GE TOL, nrej)
        IF nrej NE 0 THEN BEGIN
            max_resid = max(abs(resids[rej_inds]), jbad)
            outmask[use_inds[rej_inds[jbad]]] = 0
            use_inds = WHERE(outmask, nkeep)
        ENDIF ELSE BREAK
    ENDFOR
ENDIF ELSE xy2traceset, [0.0, 0.0], [0.0, 0.0], xset $
  , ncoeff = ncoeff, maxiter = 0 $
  , /silent, FUNC = FUNC, XMIN = XMIN, XMAX = XMAX, yfit = yfit
;; Kludge so that we get a structure out in the event of a failed fit
RETURN, xset
END
