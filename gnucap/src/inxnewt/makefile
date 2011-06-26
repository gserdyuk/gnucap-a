# makefile for newton project
include ../globalopt.mi
include ../sharedlib.mi
#include ../../thirdparty/thirdparty_tools_location.mi

OBJ_DIR	=	../objs


F_OBJS	=	\
		$(OBJ_DIR)/auxlib.o \
		$(OBJ_DIR)/linsrc.o

C_OBJS	=	\

	
CPP_OBJS =	\
		$(OBJ_DIR)/inxnslv.o \
		$(OBJ_DIR)/gmres_op.o 
	
OBJS	= 	$(F_OBJS) \
		$(C_OBJS) \
		$(CPP_OBJS)

BASE 	=	inxnewt


$(OBJ_DIR)/$(BASE).content : $(OBJS)
		echo $^ >$@

# "static rule"
$(F_OBJS): $(OBJ_DIR)/%.o: %.f		
		$(F77) -c $(GLOBALG77) $(SOFL) -o $@ $< 
		

$(C_OBJS): $(OBJ_DIR)/%.o: %.c	
		$(CC) -c $(GLOBALC) $(SOFL) -o $@ $< 
		

$(CPP_OBJS): $(OBJ_DIR)/%.o: %.cc	
		$(CPPC) -c $(GLOBALCPP) $(SOFL) -o $@ $< 

clean:
		rm $(OBJS)	
		rm $(OBJ_DIR)/$(BASE).content	

depend:
		$(MAKEDEPEND) -fmakefile -p$(OBJ_DIR)/ $(INCLUDEDIRS) *.cc
	
# DO NOT DELETE

../objs/gmres_op.o: ../API/commondefs.hh ../API/extradata.hh \
 ../API/fileio.hh ../API/mhb.hh ../API/solver.hh ../API/vector_2d.hh \
 ../include/complex_d.hh inxnewt_api.hh inxnewt_fc.hh
../objs/inxnslv.o: ../API/commondefs.hh ../API/extradata.hh \
 ../API/fileio.hh ../API/mhb.hh ../API/solver.hh ../API/vector_2d.hh \
 ../include/complex_d.hh inxnewt_api.hh inxnewt_fc.hh
