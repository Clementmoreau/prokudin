CC = c99

LFLAGS = -ltiff -ljpeg -lpng -lm

# NOTE: to change the optimization settings swap the following two lines
CFLAGS = -g
CFLAGS = -O3

all: cut3 join3 translation quantize registration main gauss pyramide irani

install: all
	cp cut3 /Users/clement/bin
	cp join3 /Users/clement/bin
	cp translation /Users/clement/bin
	cp registration /Users/clement/bin
	cp quantize /Users/clement/bin
	cp main /Users/clement/bin
	cp gauss /Users/clement/bin
	cp pyramide /Users/clement/bin
	cp irani /Users/clement/bin

cut3: cut3.c iio.c
	$(CC) $(CFLAGS) iio.c cut3.c $(LFLAGS) -o cut3

join3: join3.c iio.c
	$(CC) $(CFLAGS) iio.c join3.c $(LFLAGS) -o join3

translation: translation.c iio.c
	$(CC) $(CFLAGS) iio.c translation.c $(LFLAGS) -o translation

registration: registration.c iio.c
	$(CC) $(CFLAGS) iio.c registration.c $(LFLAGS) -o registration

quantize: quantize.c iio.c
	$(CC) $(CFLAGS) iio.c quantize.c $(LFLAGS) -o quantize

main: main.c iio.c 
	$(CC) $(CFLAGS) iio.c main.c $(LFLAGS) -o main

gauss: gauss.c iio.c
	$(CC) $(CFLAGS) iio.c gauss.c $(LFLAGS) -o gauss

pyramide: pyramide.c iio.c
	$(CC) $(CFLAGS) iio.c pyramide.c $(LFLAGS) -o pyramide

irani: irani.c iio.c
	$(CC) $(CFLAGS) iio.c irani.c $(LFLAGS) -o irani

clean:
	rm -f cut3 join3 translation quantize registration
