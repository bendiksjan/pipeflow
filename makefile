FLAGS  = -c -O5 -r8
PROF =  
LIB =  libacml.a

COMP  = /opt/openmpi/bin/mpif90 -132
#COMP  = ftn -132   aa
pijp_s:  decomp_2d.o sfft.o ders_n.o main_comp.o  mom.o solver.o output_module.o  param.txt
	$(COMP) -O5  main_comp.o mom.o decomp_2d.o io.o solver.o output_module.o sfft.o ders_n.o $(LIB) -o pijp_s
decomp_2d.o : decomp_2d.f90 
	$(COMP)  -DDOUBLE_PREC -DOVERWRITE -O3  -cpp  -c decomp_2d.f90 io.f90
main_comp.o : main_comp.f param.txt
	$(COMP) $(FLAGS) main_comp.f
mom.o : mom.f
	$(COMP) $(FLAGS) mom.f 
solver.o : solver.f param.txt
	$(COMP)  -c -132 -O5 solver.f 
sfft.o : sfft.f
	$(COMP) $(FLAGS) sfft.f
ders_n.o: ders_n.f
	$(COMP) $(FLAGS) ders_n.f
plot.o	: plof.f
	$(COMP) $(FLAGS) plot.f
output_module.o : output_module.f param.txt
	$(COMP) $(FLAGS) output_module.f
