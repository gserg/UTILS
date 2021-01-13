#          utils ION

include ${HOME}/DEFS_03

all:   bound_ovl  dif_sec_He  sdcs sdcs_omt  ddcs  tdcs   momtrans  sec_ionb

bound_ovl: bound_ovl.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR) 

tdcs: tdcs.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR) 
	rm *.mod

ddcs: ddcs.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR) 

sdcs: sdcs.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR) 

sdcs_omt: sdcs_omt.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS)  $(ZCOM)
	mv $@ $(BINDIR) 

dif_sec_He: dif_sec_He.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR) 


momtrans: momtrans.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR) 

sec_ionb: sec_ionb.f90
	$(FC) -o $@ $(link) $< $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR) 

