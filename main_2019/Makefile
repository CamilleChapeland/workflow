# Makefile for ...JMI/main

include $(CWPROOT)/src/Makefile.config

#INTEL = /opt/intel
#MKLROOT = $(INTEL)/mkl
#MKLROOT = ~/Intel2017U4/mkl
#MKLLIB=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm
MKLLIB=-mkl=sequential

L=$(CWPROOT)/lib
D = $L/libcwp.a $L/libpar.a $L/libsu.a
B = .
#B = $(DELPHIROOT)/bin
OPTC = -O2 -I$(CWPROOT)/src/Complex/include -I$(MKLROOT)/include -fPIC \
-I$(CWPROOT)/include -DMKL_ILP64 -m64 -I$(MKLROOT)/include
#OPTS NOT INCLUDED
LFLAGS= $(PRELFLAGS) -L$L -lsu -lpar -lcwp -lm $(MKLLIB) $(POSTLFLAGS)
CC=mpiicc
#g++
#-lboost_system
#LFLAGS= $(PRELFLAGS) -L$L -lsu -lpar -lcwp -lm -lboost_iostreams -lboost_filesystem $(MKLLIB) $(POSTLFLAGS)

PROGS =			\
	$B/jmi_2d_mpi \

PROGS2 = \
	jmi_2d_mpi.c subfunctions.c	\

INSTALL	:	$(PROGS)
	@-rm -f INSTALL
	@touch $@

$(PROGS):	$(CTARGET) $D 
	-$(CC) $(CFLAGS) $(PROGS2) $(LFLAGS) -o $@
	@$(MCHMODLINE)
	@echo $(@F) installed in $B

remake	:
	-rm -f $(PROGS) $(LIBS) INSTALL
	$(MAKE) 
	
clean::
	rm -f a.out junk* JUNK* core
