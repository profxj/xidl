;; Returns 1 if the input variable is undefined, 0 otherwise
FUNCTION is_undefined, var
	if size( var, /type ) EQ 0 then RETURN, 1
	RETURN,0
END
