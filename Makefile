
include make.include

all:	arg.o Grid.o StraCorn.o Skin.o LCK Calib_LCK

clean:
	rm -f *.o *.lnk

#include ../lib/make.inc
#include ../mc/make.inc


### executables ##
LCK: LCK.cpp Grid.o Skin.o
	$(LD) $(LFLAGS)  LCK.cpp Grid.o Skin.o arg.o $(LIBS) $(INC) $(LIBSPATH)/libsundials_cvode.a $(LIBSPATH)/libsundials_nvecserial.a -o LCK
Calib_LCK: Calib_LCK.cpp Grid.o Skin.o
	$(LD) $(LFLAGS)  Calib_LCK.cpp Grid.o Skin.o arg.o $(LIBS) $(INC) $(LIBSPATH)/libsundials_cvode.a $(LIBSPATH)/libsundials_nvecserial.a -o Calib_LCK

	
### objective files ##
arg.o: arg.h arg.c
	$(CC) $(CFLAGS) arg.c -o arg.o
Grid.o: Grid.h Grid.cpp
	$(CC) $(CFLAGS) Grid.cpp $(INC) -o Grid.o
StraCorn.o: StraCorn.h StraCorn.cpp
	$(CC) $(CFLAGS) StraCorn.cpp $(INC) -o StraCorn.o
Skin.o: Skin.h Skin.cpp
	$(CC) $(CFLAGS) Skin.cpp $(INC) -o Skin.o

# The following may be later used when using the MEBDFSO package
#	for large scale sparse stiff ODEs
#yale.o: yale.f
#	$(CC) $(CFLAGS) yale.f -o yale.o	
#MEBDFSO.o: MEBDFSO.f
#	$(CC) $(CFLAGS) MEBDFSO.f -o MEBDFSO.o