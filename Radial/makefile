#          utils RADIAL

include ${HOME}/DEFS_03

all:   frm_w   slater_w   w123   w_dat   w_frm

frm_w: frm_w.f90
	$(FC) -o frm_w $(link) frm_w.f90 
	mv frm_w $(BINDIR)

slater_w: slater_w.f90
	$(FC) -o slater_w $(link) slater_w.f90 $(RADIAL) $(ZCOM)
	mv slater_w $(BINDIR)

w123: w123.f90
	$(FC) -o w123 $(link) w123.f90 $(RADIAL) $(ZCOM)
	mv w123 $(BINDIR)

w_dat: w_dat.f90
	$(FC) -o w_dat $(link) w_dat.f90 $(RADIAL) $(ZCOM)
	mv w_dat $(BINDIR)

w_frm: w_frm.f90
	$(FC) -o w_frm $(link) w_frm.f90 
	mv w_frm $(BINDIR)


