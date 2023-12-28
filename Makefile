# Define compilers and flags
F90 = gfortran
CFLAGS = -g -O3
FCFLAGS = -g -O3
LDFLAGS = -L$(PWD)/libs/wigxjpf-1.11/lib
LDLIBS = -llapack -lblas -lwigxjpf -lm
INCFLAGS = -I$(PWD)/libs/wigxjpf-1.11/inc -I$(PWD)/libs/wigxjpf-1.11/mod

# Define paths for libraries
WIGXJPF_PATH = $(PWD)/libs/wigxjpf-1.11

# Define object files
OBJS = src/utility_functions_mod.o src/array_operations_mod.o src/array_operations_allocate_submod.o \
       src/array_operations_append_submod.o src/array_operations_invert_symmetric_matrix_submod.o \
       src/array_operations_fill_symmetric_matrix_submod.o src/special_functions_mod.o \
       src/math_functions_mod.o src/input_reader_mod.o src/pes_read_mod.o src/channels_mod.o \
       src/coupling_matrix_mod.o src/propagator_mod.o src/boundary_conditions_mod.o src/sts_xs_mod.o \
       src/scattering.o

.PHONY: all check_libs wigxjpf scattering

all: check_libs scattering

check_libs: wigxjpf

wigxjpf:
	@if [ ! -d "$(WIGXJPF_PATH)" ]; then \
		echo "Installing wigxjpf..."; \
		wget http://fy.chalmers.se/subatom/wigxjpf/wigxjpf-1.11.tar.gz; \
		tar -xzf wigxjpf-1.11.tar.gz; \
		cd wigxjpf-1.11; \
		make fsimple.test FC=gfortran; \
		mkdir -p $(PWD)/libs/wigxjpf-1.11/inc; \
		mkdir -p $(PWD)/libs/wigxjpf-1.11/lib; \
		mkdir -p $(PWD)/libs/wigxjpf-1.11/mod; \
		cp -r inc/* $(PWD)/libs/wigxjpf-1.11/inc/; \
		cp -r lib/* $(PWD)/libs/wigxjpf-1.11/lib/; \
		cp -r mod/* $(PWD)/libs/wigxjpf-1.11/mod/; \
		cd ..; \
		rm -rf wigxjpf-1.11 wigxjpf-1.11.tar.gz; \
	fi

scattering: $(OBJS)
	$(F90) $(CFLAGS) $(FCFLAGS) $(INCFLAGS) $^ -o scattering.x $(LDFLAGS) $(LDLIBS)
	rm src/*.o src/*.mod src/*.smod

$(OBJS) : src/%.o : src/%.f90
	$(F90) $(CFLAGS) $(FCFLAGS) $(INCFLAGS) -Jsrc -c $< -o $@

clean:
	rm -rf src/*.o src/*.mod src/*.smod *.x
