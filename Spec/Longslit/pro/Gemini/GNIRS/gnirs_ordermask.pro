FUNCTION GNIRS_ORDERMASK, newloglam, order
   ngrid = n_elements(newloglam)
   masklam = lonarr(ngrid)
   CASE order OF
       0: masklam[WHERE(10.0d^newloglam GT 19500.0D AND $
                        10.0d^newloglam LT 24000.0D)] = 1
       1: masklam[WHERE(10.0d^newloglam GT 15000.0D AND $
                        10.0d^newloglam LT 18000.0D)] = 1
       2: masklam[WHERE((10.0d^newloglam GT 12000.0D AND $
                         10.0d^newloglam LT 13500.0D) OR $
                        (10.0d^newloglam GT 14250.0  AND $
                         10.0d^newloglam  LT 15000.0))] = 1
       3: masklam[WHERE(10.0d^newloglam GT 10500.0D AND $
                        10.0d^newloglam LT 12500.0D)] = 1
       4: masklam[WHERE(10.0d^newloglam GT 9200.0D AND $
                        10.0d^newloglam LT 10800.0D)] = 1
       5: masklam[WHERE(10.0d^newloglam GT 8300.0D AND $
                        10.0d^newloglam LT 9500.0D)] = 1
   ENDCASE
   RETURN, masklam
END 
