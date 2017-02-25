
pro ps_replacetitle, replace_title, qafil

   qatmp = qafil +'.tmp'
;   spawn, "awk '{ if ($0~/%%Title/) {print " +replace_title + $
;            "} else {print $0} } ' " +  qafil + " > " + qatmp
;
;   spawn, 'mv -f ' + qatmp + ' ' + qafil

return
end
