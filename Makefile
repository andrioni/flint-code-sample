all:
	gcc -O2 -ansi -Wall -Werror -pedantic -lm qr.c -o qr

clean:
	rm qr

indent:
	indent -bap -blf -bli0 -cbi0 -cdw -cli4 -cs -i4 -l79 -nbad -nbc -nce -npcs -nprs -nut -pmt -psl -saf -sai -saw -sbi0 -ss -ts4 qr.c

debug:
	gcc -g -ansi -Wall -Werror -pedantic -lm qr.c -o qr
