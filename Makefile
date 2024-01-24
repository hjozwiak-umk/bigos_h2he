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
       src/array_operations_fill_symmetric_matrix_submod.o src/array_operations_add_scalar_to_diagonal_submod.o \
       src/special_functions_mod.o src/math_utilities_mod.o src/global_variables_mod.o src/physics_utilities_mod.o \
       src/input_validation_mod.o src/input_reader_mod.o \
       src/radial_coupling_terms_mod.o src/channels_mod.o src/pes_matrix_mod.o src/centrifugal_matrix_mod.o \
       src/propagator_mod.o src/boundary_conditions_mod.o src/unitarity_check_mod.o src/state_to_state_cross_sections_mod.o \
       src/save_s_matrix_mod.o src/scattering.o

.PHONY: all check_libs wigxjpf scattering test

all: check_libs scattering

check_libs: wigxjpf

wigxjpf:
	@if [ ! -d "$(WIGXJPF_PATH)" ]; then \
		echo "Installing wigxjpf..."; \
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
		rm -rf wigxjpf-1.11; \
	fi

scattering: $(OBJS)
	$(F90) $(CFLAGS) $(FCFLAGS) $(INCFLAGS) $^ -o scattering.x $(LDFLAGS) $(LDLIBS)
	rm src/*.o src/*.mod src/*.smod

.PHONY: test

test:
	@echo "Running test for test-elastic-oH2-He..."
	@cp scattering.x test/test-elastic-oH2-He/
	@cp data/ref/oH2-He-radialterms.zip test/test-elastic-oH2-He/
	@cd test/test-elastic-oH2-He/ && unzip oH2-He-radialterms.zip
	@cd test/test-elastic-oH2-He/ && ./scattering.x < input.dat > output.dat
	@tail -n 2 test/test-elastic-oH2-He/output.dat | head -n 1 > tmp1.txt
	@tail -n 2 data/ref/test/test-elastic-oH2-He/output.dat | head -n 1 > tmp2.txt
	@diff tmp1.txt tmp2.txt && echo "test-elastic-oH2-He passed" || echo "test-elastic-oH2-He failed"
	@rm tmp1.txt tmp2.txt test/test-elastic-oH2-He/oH2-He-radialterms* test/test-elastic-oH2-He/scattering.x

	@echo "Running test for test-orthoH2-He..."
	@cp scattering.x test/test-orthoH2-He/
	@cp data/ref/oH2-He-radialterms.zip test/test-orthoH2-He/
	@cd test/test-orthoH2-He/ && unzip oH2-He-radialterms.zip
	@cd test/test-orthoH2-He/ && ./scattering.x < input.dat > output.dat
	@tail -n 37 test/test-orthoH2-He/output.dat | head -n 36 > tmp1.txt
	@tail -n 37 data/ref/test/test-orthoH2-He/output.dat | head -n 36 > tmp2.txt
	@diff tmp1.txt tmp2.txt && echo "test-orthoH2-He passed" || echo "test-orthoH2-He failed"
	@rm tmp1.txt tmp2.txt test/test-orthoH2-He/oH2-He-radialterms* test/test-orthoH2-He/scattering.x
	
	@echo "Running test for test-paraH2-He..."
	@cp scattering.x test/test-paraH2-He/
	@cp data/ref/pH2-He-radialterms.zip test/test-paraH2-He/
	@cd test/test-paraH2-He/ && unzip pH2-He-radialterms.zip
	@cd test/test-paraH2-He/ && ./scattering.x < input.dat > output.dat
	@tail -n 226 test/test-paraH2-He/output.dat | head -n 225 > tmp1.txt
	@tail -n 226 data/ref/test/test-paraH2-He/output.dat | head -n 225 > tmp2.txt
	@diff tmp1.txt tmp2.txt && echo "test-paraH2-He passed" || echo "test-paraH2-He failed"
	@rm tmp1.txt tmp2.txt test/test-paraH2-He/pH2-He-radialterms* test/test-paraH2-He/scattering.x

$(OBJS) : src/%.o : src/%.f90
	$(F90) $(CFLAGS) $(FCFLAGS) $(INCFLAGS) -Jsrc -c $< -o $@

clean:
	rm -rf src/*.o src/*.mod src/*.smod *.x
