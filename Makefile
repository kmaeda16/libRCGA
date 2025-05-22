CC = /usr/bin/gcc # CHANGE HERE
MPICC = /usr/bin/mpicc # CHANGE HERE
CFLAGS = -Wall -O2

INCDIR = include
SRCDIR = source

TGT1 = libundxmgg_serial
RCGATYPE1 = UNDX_MGG_serial

TGT2 = librexstarjgg_serial
RCGATYPE2 = REXstar_JGG_serial

TGT3 = libundxmgg_parallel
RCGATYPE3 = UNDX_MGG_parallel

TGT4 = librexstarjgg_parallel
RCGATYPE4 = REXstar_JGG_parallel


all: lib/$(TGT1).a lib/$(TGT2).a lib/$(TGT3).a lib/$(TGT4).a


# ---------------- libundxmgg_serial.a ---------------- #

OBJS1 = $(SRCDIR)/$(RCGATYPE1).o $(SRCDIR)/Shared.o
PREREQS1 = $(SRCDIR)/$(RCGATYPE1).c $(SRCDIR)/Shared.c \
           $(INCDIR)/$(RCGATYPE1).h $(INCDIR)/Shared.h

lib/$(TGT1).a: $(OBJS1)
	(ar rcv lib/$(TGT1).a $(OBJS1); ranlib lib/$(TGT1).a)

$(SRCDIR)/$(RCGATYPE1).o: $(PREREQS1)
	(cd $(SRCDIR); $(CC) $(CFLAGS) -I../$(INCDIR) -c $(RCGATYPE1).c)


# ---------------- librexstarjgg_serial.a ---------------- #

OBJS2 = $(SRCDIR)/$(RCGATYPE2).o $(SRCDIR)/Shared.o
PREREQS2 = $(SRCDIR)/$(RCGATYPE2).c $(SRCDIR)/Shared.c \
           $(INCDIR)/$(RCGATYPE2).h $(INCDIR)/Shared.h

lib/$(TGT2).a: $(OBJS2)
	(ar rcv lib/$(TGT2).a $(OBJS2); ranlib lib/$(TGT2).a)

$(SRCDIR)/$(RCGATYPE2).o: $(PREREQS2)
	(cd $(SRCDIR); $(CC) $(CFLAGS) -I../$(INCDIR) -c $(RCGATYPE2).c)


# ---------------- libundxmgg_parallel.a ---------------- #

OBJS3 = $(SRCDIR)/$(RCGATYPE3).o $(SRCDIR)/Shared.o
PREREQS3 = $(SRCDIR)/$(RCGATYPE3).c $(SRCDIR)/Shared.c \
           $(INCDIR)/$(RCGATYPE3).h $(INCDIR)/Shared.h

lib/$(TGT3).a: $(OBJS3)
	(ar rcv lib/$(TGT3).a $(OBJS3); ranlib lib/$(TGT3).a)

$(SRCDIR)/$(RCGATYPE3).o: $(PREREQS3)
	(cd $(SRCDIR); $(MPICC) $(CFLAGS) -I../$(INCDIR) -c $(RCGATYPE3).c)


# ---------------- librexstarjgg_parallel.a ---------------- #

OBJS4 = $(SRCDIR)/$(RCGATYPE4).o $(SRCDIR)/Shared.o
PREREQS4 = $(SRCDIR)/$(RCGATYPE4).c $(SRCDIR)/Shared.c \
           $(INCDIR)/$(RCGATYPE4).h $(INCDIR)/Shared.h

lib/$(TGT4).a: $(OBJS4)
	(ar rcv lib/$(TGT4).a $(OBJS4); ranlib lib/$(TGT4).a)

$(SRCDIR)/$(RCGATYPE4).o: $(PREREQS4)
	(cd $(SRCDIR); $(MPICC) $(CFLAGS) -I../$(INCDIR) -c $(RCGATYPE4).c)


# ---------------- Shared.o ---------------- #

$(SRCDIR)/Shared.o: $(SRCDIR)/Shared.c $(INCDIR)/Shared.h
	(cd $(SRCDIR); $(CC) $(CFLAGS) -I../$(INCDIR) -c Shared.c)