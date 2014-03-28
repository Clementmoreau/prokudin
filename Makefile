CC = c99

LFLAGS = -ltiff -ljpeg -lpng -lm

# NOTE: to change the optimization settings swap the following two lines
CFLAGS = -O3
CFLAGS = -g

ALL = cut3 join3 translation quantize registration main gauss pyramide irani echantillon

all: $(ALL)

install: all
	cp $(ALL) /Users/clement/bin

$(ALL) : % : %.c iio.o
	$(CC) $(CFLAGS) $^ $(LFLAGS) -o $@

clean:
	rm -f $(ALL) iio.o
