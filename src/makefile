F90 = gfortran # compile with gfortran
########################################################################
# wigxjpf library 
########################################################################
# H. T. Johansson and C. Forssén,
# Fast and Accurate Evaluation of Wigner 3j, 6j, and 9j Symbols Using 
# Prime Factorization and Multiword Integer Arithmetic, 
# SIAM J. Sci. Comput., 38(1) (2016), A376-A384. 
# eprint http://dx.doi.org/10.1137/15M1021908
# Pre-print (2015) at arXiv:1504.08329
########################################################################
wig_prefix=/home/hubert/Work/H2-He/RateCoeffs/2022/paper/review/bigos/libs/wigxjpf-1.11
wig_libdir=${wig_prefix}/lib
wig_includedir=${wig_prefix}/inc
wig_moduledir=${wig_prefix}/mod
wig_name=wigxjpf
########################################################################
CMPLFLG = -c -O3 -I${wig_includedir} -I${wig_moduledir} 
Libs = -L${wig_libdir} -l${wig_name}
########################################################################

ADD_OBJS = utility_functions_mod.o array_operations_mod.o array_operations_allocate_submod.o array_operations_append_submod.o array_operations_invert_symmetric_matrix_submod.o array_operations_fill_symmetric_matrix_submod.o
OBJS = $(ADD_OBJS)  special_functions_mod.o math_functions_mod.o input_reader_mod.o pes_read_mod.o channels_mod.o coupling_matrix_mod.o propagator_mod.o boundary_conditions_mod.o sts_xs_mod.o scattering.o
.SUFFIXES: .f .f90 .o

all : $(OBJS)
	$(F90) $(OBJS) -o scattering.x $(Libs) -llapack -lblas -Wunused-variable
	rm *.o 
	rm *.mod 

clean:
	rm -rf *.o *.mod *.smod *.x

$(OBJS) :
.f90.o:
	$(F90) $(CMPLFLG) $<
.f.o:
	$(F90) $(CMPLFLG) $<
