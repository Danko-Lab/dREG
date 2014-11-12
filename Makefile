R_LIBS := ~/bin/dREG
export R_LIBS

R_dependencies:
	@echo "Installing R dependencies to:" ${R_LIBS}
	mkdir -p ${R_LIBS}
	R --no-save < rDeps.R

dreg:
	@echo "Installing dREG to:" ${R_LIBS}
	mkdir -p ${R_LIBS}
	make topLibs -C kent/src
	R CMD INSTALL bigWig --clean
	make -C dREG/src/lib
	R CMD INSTALL dREG --clean

uninstall: 
	rm -Rf ${R_LIBS}

