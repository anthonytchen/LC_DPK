# C-compiler
CC=g++

# libraries
LIBS=-L/usr/local/lib -lm -pthread -lsundials_cvode -lsundials_nvecserial
INC=-I/usr/local/include/

#CFLAGS=-O2 -ffast-math -Wno-deprecated -c
CFLAGS=-g -O0 -c -Wno-write-strings


# linker
LD=g++

#LFLAGS=-O2 -ffast-math -Wno-deprecated
LFLAGS=-g -O0 -Wno-write-strings


all:	except.o arg.o Chemical.o \
	Grid.o Comp.o Vehicle.o Sebum.o StraCorn.o ViaEpd.o Dermis.o Skin.o Skin_VS.o Skin_VSVDB.o Blood.o \
	run_VS run_VSVDB #LCK_metabo

clean:
	rm -f *.o *.lnk

### executables ##

# vehicle & stratum corneum
run_VS: run_VS.cpp except.o arg.o Chemical.o \
	Grid.o Comp.o Vehicle.o Sebum.o StraCorn.o ViaEpd.o Dermis.o Skin.o Skin_VS.o Blood.o
	$(LD) $(LFLAGS)  run_VS.cpp except.o arg.o Chemical.o \
	Grid.o Comp.o Vehicle.o Sebum.o StraCorn.o ViaEpd.o Dermis.o Skin.o Skin_VS.o Blood.o \
	$(LIBS) $(INC) -o run_VS

# vehicle, stratum corneum, viable epidermis, dermis & blood
run_VSVDB: run_VSVDB.cpp except.o arg.o Chemical.o \
	Grid.o Comp.o Vehicle.o Sebum.o StraCorn.o ViaEpd.o Dermis.o Skin.o Skin_VSVDB.o Blood.o
	$(LD) $(LFLAGS)  run_VSVDB.cpp except.o arg.o Chemical.o \
	Grid.o Comp.o Vehicle.o Sebum.o StraCorn.o ViaEpd.o Dermis.o Skin.o Skin_VSVDB.o Blood.o \
	$(LIBS) $(INC) -o run_VSVDB

#LCK_metabo: LCK_metabo.cpp Chemical.o Grid.o StraCorn.o ViaEpd.o Dermis.o Skin.o Blood.o
#	$(LD) $(LFLAGS)  LCK_metabo.cpp Chemical.o Grid.o StraCorn.o ViaEpd.o Dermis.o Skin.o Blood.o arg.o $(LIBS) $(INC) -o LCK_metabo



### objective files ##

except.o: except.cpp except.h
	$(CC) $(CFLAGS) except.cpp $(INC) -o except.o
arg.o: arg.h arg.c
	$(CC) $(CFLAGS) arg.c $(INC) -o arg.o
Chemical.o: Chemical.h Chemical.cpp
	$(CC) $(CFLAGS) Chemical.cpp $(INC) -o Chemical.o
Grid.o: Grid.h Grid.cpp
	$(CC) $(CFLAGS) Grid.cpp $(INC) -o Grid.o
Comp.o: Comp.h Comp.cpp
	$(CC) $(CFLAGS) Comp.cpp $(INC) -o Comp.o
Vehicle.o: Vehicle.h Vehicle.cpp
	$(CC) $(CFLAGS) Vehicle.cpp $(INC) -o Vehicle.o
Sebum.o: Sebum.h Sebum.cpp
	$(CC) $(CFLAGS) Sebum.cpp $(INC) -o Sebum.o
StraCorn.o: StraCorn.h StraCorn.cpp
	$(CC) $(CFLAGS) StraCorn.cpp $(INC) -o StraCorn.o
ViaEpd.o: ViaEpd.h ViaEpd.cpp
	$(CC) $(CFLAGS) ViaEpd.cpp $(INC) -o ViaEpd.o
Dermis.o: Dermis.h Dermis.cpp
	$(CC) $(CFLAGS) Dermis.cpp $(INC) -o Dermis.o
Skin.o: Skin.h Skin.cpp
	$(CC) $(CFLAGS) Skin.cpp $(INC) -o Skin.o
Skin_VS.o: Skin_VS.h Skin_VS.cpp
	$(CC) $(CFLAGS) Skin_VS.cpp $(INC) -o Skin_VS.o
Skin_VSVDB.o: Skin_VSVDB.h Skin_VSVDB.cpp
	$(CC) $(CFLAGS) Skin_VSVDB.cpp $(INC) -o Skin_VSVDB.o
Blood.o: Blood.h Blood.cpp
	$(CC) $(CFLAGS) Blood.cpp $(INC) -o Blood.o
