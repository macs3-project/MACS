CC=			gcc
CFLAGS=		-g -Wall -O2 -Wno-unused-function #-fno-inline-functions -fno-inline-functions-called-once
CPPFLAGS=
INCLUDES=	
OBJS=		kthread.o misc.o \
			bseq.o htab.o bfc.o \
			rle.o rope.o mrope.o rld0.o \
			unitig.o mag.o bubble.o ksw.o
PROG=		fml-asm
LIBS=		-lm -lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

fml-asm:libfml.a example.o
		$(CC) $(CFLAGS) $^ -o $@ -L. -lfml $(LIBS)

libfml.a:$(OBJS)
		$(AR) -csru $@ $(OBJS)

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c)

# DO NOT DELETE

bfc.o: htab.h kmer.h internal.h fml.h kvec.h ksort.h
bseq.o: fml.h kseq.h
bubble.o: mag.h kstring.h fml.h kvec.h ksw.h internal.h khash.h
example.o: fml.h
htab.o: htab.h kmer.h khash.h
ksw.o: ksw.h
mag.o: mag.h kstring.h fml.h kvec.h internal.h kseq.h khash.h ksort.h
misc.o: internal.h fml.h kstring.h rle.h mrope.h rope.h rld0.h mag.h kvec.h
misc.o: htab.h kmer.h khash.h
mrope.o: mrope.h rope.h
rld0.o: rld0.h
rle.o: rle.h
rope.o: rle.h rope.h
unitig.o: kvec.h kstring.h rld0.h mag.h fml.h internal.h ksort.h
