########################################################################
#                  Makefile to compile NIKA
IDIR =./src                                        #paramters, constants
BDIR =../bin                        #will be the final binary executable
ODIR = obj                                        #subroutines in binary 
SDIR =.                                                     #orio akì xd
LIBS = -llapack                                            #library name
gg = gfortran -I$(IDIR)                                        #compiler

ifeq ($(DEBUG),1)
FFLAGS = -Wall -g -msse4.2 -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant
else
FFLAGS = -Wall -Wno-unused -Wno-unused-dummy-argument -O2
endif

SRCF90 = $(wildcard *.f90)
OBJ = $(patsubst %.f90,$(ODIR)/%.o,$(SRCF90))

$(ODIR)/%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS) $(LIBS)

$(BDIR)/hfmp2dft: $(OBJ)
	$(gg) -o $@.exe $^ $(FFLAGS) $(LIBS)
	@echo " ----------------------------------------------------- "
	@echo "   ✨🌟 HFMP2DFT has been successfully compiled 🌟✨   "
#	@echo " "
#	@echo "             run: ./launcher.sh input"
	@echo " ----------------------------------------------------- "

debug:
	DEBUG=1 make clean $(BDIR)/hfmp2dft
clean:
	rm -f $(ODIR)/*.o $(BDIR)/hfmp2dft $(BDIR)/debug
#	rm -f $(OBJ).exe *.o *.mod *.a

