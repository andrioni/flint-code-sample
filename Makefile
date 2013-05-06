CC=gcc
CFLAGS=-O3 -ansi -Wall -Werror -pedantic
DEBUGFLAGS=-g -ansi -Wall -Werror -pedantic
LDFLAGS=-lm -lmpfr
FILES=util.c qr.c
COLFILES=qr_col.c
HEADERS=qr.h

all: normal col

normal:
	$(CC) $(CFLAGS) $(LDFLAGS) $(FILES) -o qr

col:
	$(CC) $(CFLAGS) $(LDFLAGS) $(COLFILES) -o qr_col

noprint:
	$(CC) -DPRINTOUTPUT=0 $(CFLAGS) $(LDFLAGS) $(FILES) -o qr
	$(CC) -DPRINTOUTPUT=0 $(CFLAGS) $(LDFLAGS) $(COLFILES) -o qr_col

clean:
	rm qr
	rm qr_col
	rm .c~

indent:
	indent -bap -blf -bli0 -cbi0 -cdw -cli4 -cs -i4 -l79 -nbad -nbc -nce -npcs -nprs -nut -pmt -psl -saf -sai -saw -sbi0 -ss -ts4 $(FILES) $(COLFILES) $(HEADERS)

debug:
	$(CC) $(DEBUGFLAGS) $(LDFLAGS) $(FILES) -o qr

debugcol:
	$(CC) $(DEBUGFLAGS) $(LDFLAGS) $(COLFILES) -o qr

bench: noprint
	time ./qr < test/100x100.data
	time ./qr < test/500x500.data
	time ./qr < test/1000x1000.data
	time ./qr_col < test/100x100.data
	time ./qr_col < test/500x500.data
	time ./qr_col < test/1000x1000.data
