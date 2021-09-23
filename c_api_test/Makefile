CFLAGS += -g -I.. -I../htslib/htslib
# Required to link with g++ because libxsqueezeit is compiled with g++
LD = g++

LDFLAGS = -g -L ../htslib
LIBS=-lpthread -lhts -lzstd

all : main

main : main.o libxsqueezeit.a
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

libxsqueezeit.a : xsqueezeit_export
	cp ../xsqueezeit_export/libxsqueezeit.a .

xsqueezeit_export :
	cd ..; make package-sources
	cd ../xsqueezeit_export; make

clean :
	rm -rf *.o
	rm -rf *.d
	rm -f main
	rm -f libxsqueezeit.a

.PHONY :
	all clean xsqueezeit_export