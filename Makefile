
include make.include

all:	arg.o Chemical.o Grid.o StraCorn.o ViaEpd.o Dermis.o Skin.o Blood.o LCK

clean:
	rm -f *.o *.lnk

#include ../lib/make.inc
#include ../mc/make.inc


### executables ##
LCK: LCK.cpp Chemical.o Grid.o StraCorn.o ViaEpd.o Dermis.o Skin.o Blood.o
	$(LD) $(LFLAGS) LCK.cpp Chemical.o Grid.o StraCorn.o ViaEpd.o Dermis.o Skin.o Blood.o arg.o $(LIBS) $(INC) $(LIBSPATH)/libsundials_cvode.a $(LIBSPATH)/libsundials_nvecserial.a -o LCK
Calib_LCK: Calib_LCK.cpp Grid.o Skin.o
	$(LD) $(LFLAGS)  Calib_LCK.cpp Grid.o Skin.o arg.o $(LIBS) $(INC) $(LIBSPATH)/libsundials_cvode.a $(LIBSPATH)/libsundials_nvecserial.a -o Calib_LCK


### objective files ##
arg.o: arg.h arg.c
	$(CC) $(CFLAGS) arg.c -o arg.o
Chemical.o: Chemical.h Chemical.cpp
	$(CC) $(CFLAGS) Chemical.cpp $(INC) -o Chemical.o
Grid.o: Grid.h Grid.cpp
	$(CC) $(CFLAGS) Grid.cpp $(INC) -o Grid.o
StraCorn.o: StraCorn.h StraCorn.cpp
	$(CC) $(CFLAGS) StraCorn.cpp $(INC) -o StraCorn.o
ViaEpd.o: ViaEpd.h ViaEpd.cpp
	$(CC) $(CFLAGS) ViaEpd.cpp $(INC) -o ViaEpd.o
Dermis.o: Dermis.h Dermis.cpp
	$(CC) $(CFLAGS) Dermis.cpp $(INC) -o Dermis.o
Skin.o: Skin.h Skin.cpp
	$(CC) $(CFLAGS) Skin.cpp $(INC) -o Skin.o
Blood.o: Blood.h Blood.cpp
	$(CC) $(CFLAGS) Blood.cpp $(INC) -o Blood.o
