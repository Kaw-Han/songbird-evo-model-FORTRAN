# This defines the name of the main executable file
EXEC_NAME=model

# Defines the file name for the test program
TESTS_NAME=tests

# Output data from simulation: data for each generation
OUTPUT_CSV=my_model_output.csv

# The list of plot files that are generated from the output data
PLOTS=plot_gens-1.svg plot_gens-2.svg plot_gens-3.svg plot_gens-4.svg plot_gens-5.svg plot_gens-6.svg

# These makefile macros define which compiler flags are
# used to build the target(s). The FFLAGS (abbreviation of Fortran FLAGS)
# variable keeps the parameters.
#
# If the variable DEBUG is defined to any value, then the fortran compiler
# will insert debug symbols into the code (-g) and switch off any compiler
# optimizations (-O0) - this makes it easier to locate errors
ifdef DEBUG
	FFLAGS = -g -O0
endif

# If the DEBUG variable is not defines, the fortran compiler uses the highest
# optimization (-O3) and fast mathermatical optimization (-ffast-math), so the
# resulting executable program (usually) runs much faster.
ifndef DEBUG
  FFLAGS = -O3 -ffast-math
endif

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Determine what is the build platform, Windows / non-Windows
#
# A safer way to check platform if uname is not available, ComSpec on Windows
# Note that ComSpec is (may be?) case-sensitive, check with env.exe
ifdef ComSpec
	PLATFORM_TYPE=Windows
	WHICH_CMD=where
	NULLDEV=":NULL"
	ECHO := "$(shell $(WHICH_CMD) echo.exe)"
	RM := del /Q
	MV := move
else
	PLATFORM_TYPE=Unix
	WHICH_CMD=which
	NULLDEV="/dev/null"
	ECHO := echo
	RM := rm -fr
	RMDIR := rmdir
	MV := mv -f
endif

# --- TARGETS FOLLOW BELOW -----

# The all target is the main default target that is called
# with 'make' or 'make all':
# It should refer (after : sign) to what is built by default, e.g.
#    all: $(TESTS_NAME) to make the tests
#    all: $(EXEC_NAME) to make the main simulation program
#    all: $(OUTPUT_CSV) to run the simulation and produce the main output file
#    all: $(PLOTS) to run the simulation and produce the plots
all: $(PLOTS)

# Run produces the output CSV file
run: $(OUTPUT_CSV)

plots: $(PLOTS)

$(PLOTS): $(OUTPUT_CSV)
	gnuplot plot.gnuplot

# This target builds the main model executable from the source code
$(EXEC_NAME): m_utils.f90 m_random.f90 m_csv_data.f90 first_model.f90 main.f90
	gfortran $(FFLAGS) -o $(EXEC_NAME) $^

# This target builds the test program executable from the source code
$(TESTS_NAME): m_utils.f90 m_random.f90 m_csv_data.f90 first_model.f90 tests.f90
	gfortran -g -O0 -o $(TESTS_NAME) $^

	# Run the main simulation
$(OUTPUT_CSV): $(EXEC_NAME)
	./$(EXEC_NAME)

# Cleanup all: delete all temporary files and programs
clean: cleandata
	$(RM) *.mod $(EXEC_NAME) $(EXEC_NAME).exe $(TESTS_NAME) $(TESTS_NAME).exe *.tmp

cleandata:
	-$(RM) generations/* $(OUTPUT_CSV) $(PLOTS)
	-$(RMDIR) generations
