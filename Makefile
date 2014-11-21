RSLPATH = /usr/local/trmm
INCDIR = -I$(RSLPATH)/include
LIBFLAGS =  -L/usr/lib -L$(RSLPATH)/lib -Wl,-rpath=$(RSLPATH)/lib -lm 

CC = gcc

CFLAGS =  -g -I. $(INCDIR) -I.. 

all: 
	$(CC) -fPIC rsl_to_hdf5odim.c -c  $(CFLAGS) $(LIBFLAGS)
	$(CC)  -shared -o librsl_to_hdf5odim.so  rsl_to_hdf5odim.o
