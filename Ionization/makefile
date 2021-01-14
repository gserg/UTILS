#          utils ION

include ${HOME}/DEFS_03

all:   bound_ovl  ddcs  dif_sec_He  momtrans sdcs sdcs_omt sec_ionb tdcs     

bound_ovl: bound_ovl.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR) 

ddcs: ddcs.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR) 

dif_sec_He: dif_sec_He.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR) 

momtrans: momtrans.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR) 

sdcs: sdcs.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR) 

sdcs_omt: sdcs_omt.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS)  $(ZCOM)
	mv $@ $(BINDIR) 

tdcs: tdcs.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR) 
	rm *.mod


sec_ionb: sec_ionb.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR) 

