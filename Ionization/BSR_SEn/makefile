include ${HOME}/DEFS_03
 
N  = bsr_sen1
 
# Libraries to be searched:
 
L  =   $(ZCONFLS) $(ZCOM) $(LAPACK)

S  =	mod_ovl.f90            bsr_sen.f90         read_data.f90   \
	read_Hdat.f90          rm_ovlF.f90         rovl_out.f90    \
        sub_sct.f90            zaface.f90  
 
O = $(S:.f90=.o) crees.o

 
$(N): $(O)
	$(FC) -o $(N) $(O) $(link) $(L)
	mv $(N) $(BINDIR)
	rm *.o *.mod 
 
clean:
	rm -f *.o *.mod 
 
x:
	rm -f $(N)
 
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
 
.f90.o:
	$(FC) $(comp)  $<

crees.o: crees.f
	$(FC) $(comp) crees.f
