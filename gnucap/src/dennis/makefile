include ../globalopt.mi
include ../sharedlib.mi

BASE	=	dennis
LIBNAME	=	$(BASE).a
OBJ_DIR	=	../objs


DENNIS_SRC	= \
		newsol.f \
		stop.f \
		stop0.f \
		mache.f \
		gradie.f \
		lineq1.f \
		lsearch.f \
		neinck.f \
		fdjac1_d.f \
		enorm.f \
		dpmpar.f 
		
DENNIS_OBJS	= \
		$(OBJ_DIR)/newsol.o \
		$(OBJ_DIR)/stop.o \
		$(OBJ_DIR)/stop0.o \
		$(OBJ_DIR)/mache.o \
		$(OBJ_DIR)/gradie.o \
		$(OBJ_DIR)/lineq1.o \
		$(OBJ_DIR)/lsearch.o \
		$(OBJ_DIR)/neinck.o \
		$(OBJ_DIR)/fdjac1_d.o \
		$(OBJ_DIR)/enorm.o \
		$(OBJ_DIR)/dpmpar.o
		

F_OBJS	=	\
		$(DENNIS_OBJS)

C_OBJS	=	\
	
CPP_OBJS =	\

OBJS	= 	$(F_OBJS)\
		$(C_OBJS) \
		$(CPP_OBJS)

FFLAGS	= $(GLOBALG77) 	$(SOFL)
CFLAGS 	= $(GLOBALC) $(SOFL)
CPPFLAGS= $(GLOBALCPP) $(SOFL)

$(OBJ_DIR)/$(BASE).content:	$(OBJS)	
	echo $(OBJS) >$(OBJ_DIR)/$(BASE).content
		
$(LIBNAME) :  $(OBJS)
	$(LIBR) -rc $@ $(OBJS)
	$(RANLIB) $@


# "static rule"
$(F_OBJS): $(OBJ_DIR)/%.o: %.f		
		$(F77) -c $(FFLAGS) -o $@ $< 
		

$(C_OBJS): $(OBJ_DIR)/%.o: %.c	
		$(CC) -c $(CFLAGS)  -o $@ $< 
		

$(CPP_OBJS): $(OBJ_DIR)/%.o: %.c	
		$(CPPC) -c $(CPPFLAGS)  -o $@ $< 

clean:
	rm -f $(OBJ_DIR)/$(BASE).content
	rm -f $(OBJS)	
	

depend:
	
