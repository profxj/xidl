pro parsehfiles

; This algorithm searches through all the .h files for definitions of 
; system variables, then writes an IDL routine called initdeimos.pro that
; sets up those variables

spawn,"grep -h define *.h | grep -v #define | grep -v '# define' | grep -v '$1' > grepdef.txt"

openr,2,"grepdef.txt"
openw,3,"initdeimos.pro"
printf,3,"pro initdeimos"
printf,3,";This routine initializes system variables for the optical model."
printf,3
printf,3

while not eof(2) do begin
	tmp=""
	readf,2,tmp

	tmp=strcompress(strtrim(tmp,2))
	tmps=strsplit(tmp,/extract)

	printf,3,"defsysv,'!",tmps(1),"',",tmps(2)
endwhile

close,2

printf,3
printf,3

printf,3,"return"
printf,3,"end"
close,3



return
end
