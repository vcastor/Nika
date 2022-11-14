########################################################################
#              Nika project Make file,
#                           to just type make 🚀
IDIR =../include
BDIR =../bin
ODIR = obj
SDIR =.
gg = gfortran -I$(IDIR)                                        #compiler
LIBS = -llapack -lblas                                        #libreries
FLAG = -fcheck=all -Wall -O3

#common = $(wildcard common_sub/*.f90)
#dft = $(wildcard dft/*.f90)
#hf = $(wildcard hf/*.f90)
#integrals = $(wildcard integrals/*.f90)
#mp2 = $(wildcard mp2/*.f90)

#objcommon = $(patsubst %.f90,$(ODIR)/%.o,$(common))
#OBJ = $(objcommon) $(objdft) $(objhf) $(objintegrals) $(objmp2)

#objcommon = $(patsubst %.f90,$(ODIR)/%.o,$(common))

SRCF90 = $(wildcard ./*/*.f90) $(wildcard *.f90)

OBJ = $(patsubst %.f90,$(ODIR)/%.o,$(SRCF90))


$(ODIR)/%.o: %.f90
	$(gg) -c -o $@ $< $(FLAG)

$(BDIR)/nika.exe: $(OBJ)
	$(gg) -o $@ $^ $(LIBS) $(FLAG)
	@echo " ------------------------------------------------------ "
	@echo "         ✨ Nika has been compiled correctly ✨"
	@echo " ------------------------------------------------------ "


clean:
	rm -f $(ODIR)/*.o $(ODIR)/*/*.o $(BDIR)/*.exe
