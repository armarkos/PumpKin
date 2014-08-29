OBJ = PumpKin.o pk_IO.o pk_Array.o pk_Pathways.o pk_Print.o
#Auxiliary.o Tests.o Print.o Pathways.o Optimize.o OptimizeDifferent.o

#CFLAGS    = -std=gnu++0x -g 
DFLAGS    = -DCOMPILING_HOSTNAME=\"$(HOSTNAME)\"
INCLUDE_DIRS += -I/opt/local/include -I$(HOME)/include
LFLAGS    =  -stdlib=libc++
LIBS =  -L/opt/local/lib -lglpk -lm 
#LIBS += = -L/usr/local/NLopt/lib -lnlopt -lm
CC   = g++
LINK = g++

all:	pumpkin

%.o:	%.cpp
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDE_DIRS) -O3 -o $@ -c $< 
#	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDE_DIRS) -g -o $@ -c $< 

%.o:	%.cc
	$(CC) $(CFLAGS) $(DFLAGS) $(INCLUDE_DIRS) -g -o $@ -c $< 

pumpkin:  $(OBJ)
	$(LINK) $(LFLAGS) -std=gnu++0x $(OBJ) -o $@ $(LIBS)


clean:  
	rm *.o *.*~ pumpkin 2> /dev/null

