# Compiler
FC = pgfortran

#Flags
FCFLAGS = -O3
HOMEPATH=$(PWD)

grem: fortran/grem_sort.f90 fortran/grem_analyze.f90
	@echo "Beginning gREM build"
	@echo ${HOMEPATH}
	mkdir -p bin/
	touch bin/test
	rm bin/*
	$(FC) $(FCFLAGS) -o bin/grem_sort.exe fortran/grem_sort.f90
	$(FC) $(FCFLAGS) -o bin/grem_analyze.exe fortran/grem_analyze.f90
	ln -s $(HOMEPATH)/python/auto_gen.py bin/
	ln -s $(HOMEPATH)/python/dump_analyze.py  bin/
	ln -s $(HOMEPATH)/python/maxwell-construction.py  bin/
	ln -s $(HOMEPATH)/python/pull_configs.py  bin/
	ln -s $(HOMEPATH)/python/read_log.py  bin/
	ln -s $(HOMEPATH)/python/run_gREM.py  bin/
	ln -s $(HOMEPATH)/python/set_analysis.py  bin/
	ln -s $(HOMEPATH)/python/sortdumps.py bin/
	chmod 777 bin/*



