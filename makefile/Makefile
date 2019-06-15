#MAIN = ../example/LaplaceEigen
#MAIN = ../example/LaplaceEigenP2

#MAIN = ../example/simplest







SELFLIB = ../lib/solvers.a
include ./make.inc
# If you want to check the momory, please use 
# valgrind --leak-check=full ./main > check

$(MAIN) : $(OBJECTS)
	$(LINKER)  $(CPPFLAGS) $(LDFLAGS)  $(CFLAGS) -fopenmp  -o $@ $^ $(LIBS)
.c.o:
	$(CC) $(CPPFLAGS)  $(CFLAGS)  $(LIBS) -fopenmp  -o $*.o -c $*.c


.PHONY:lib
lib : $(SELFLIB)
$(SELFLIB) : $(filter-out $(MAIN).o, $(OBJECTS))
	ar cr $@ $^

	
clean : 
	rm -f $(OBJECTS) 
	rm -rf $(SELFLIB)
	rm -rf $(MAIN)
	rm -rf *~
	rm -rf ../*~
	rm -rf ../*/*~
	rm -rf ../*/*/*~
