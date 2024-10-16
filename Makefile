# This defines the name of the main executable file
EXEC_NAME=model

# Defines the file name for the test program
TESTS_NAME=tests

# Defines the file name for the emotion_step test program
EMOTION_STEP_NAME=test_emotion

# Output data from simulation: data for each generation
OUTPUT_CSV = my_model_output_2G.csv
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

# Define simulation folders
SIM_FOLDERS = evol_sim_1 evol_sim_2 evol_sim_3 evol_sim_4 evol_sim_5 evol_sim_6 evol_sim_7 evol_sim_8 evol_sim_9 evol_sim_10
SIM_FOLDERS_TEST = evol_test_1
# The all target is the main default target that is called
all: run_test

# Run 1 test simulation "make" in terminal
run_test: $(EXEC_NAME)
	$(ECHO) "Running test simulation in current directory"
	./$(EXEC_NAME)
# Run 10 simulations
run_simulations: $(EXEC_NAME)
	@for folder in $(SIM_FOLDERS); do \
		mkdir -p $$folder; \
		$(ECHO) "Running simulation in $$folder"; \
		./$(EXEC_NAME) > $$folder/simulation_output.log; \
		$(MV) $(OUTPUT_CSV) $$folder/; \
		$(MV) generations $$folder/; \
		$(MV) weight_history_2G $$folder/; \
		$(MV) emotion_history_2G $$folder/; \
		$(MV) gene_history_2G $$folder/; \
	done
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

# This target builds the emotion_step testing program executable from the source code
$(EMOTION_STEP_NAME): m_utils.f90 m_random.f90 m_csv_data.f90 first_model.f90 test_emotion_step.f90
	gfortran -g -O0 -o $(EMOTION_STEP_NAME) $^
	./$(EMOTION_STEP_NAME)

# Run the main simulation
$(OUTPUT_CSV): $(EXEC_NAME)
	./$(EXEC_NAME)

# Cleanup all: delete all temporary files and programs
clean: cleandata
	$(RM) *.mod $(EXEC_NAME) $(EXEC_NAME).exe $(TESTS_NAME) $(EMOTION_STEP_NAME) $(TESTS_NAME).exe *.tmp *.csv
	@for folder in $(SIM_FOLDERS); do \
		$(RM) -r $$folder; \
	done

cleandata:
	-$(RM) generations/* $(OUTPUT_CSV) $(PLOTS)
	-$(RMDIR) generations
	-$(RM) weight_history_2G/* $(OUTPUT_CSV)
	-$(RMDIR) weight_history_2G
	-$(RM) emotion_history_2G/* $(OUTPUT_CSV)
	-$(RMDIR) emotion_history_2G
	-$(RM) gene_history_2G/* $(OUTPUT_CSV)
	-$(RMDIR) gene_history_2G/
