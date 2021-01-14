#       utilities DBSR

include ${HOME}/DEFS_03

all: bsw_rw       dbound_bsw    dbound_tab    dbsr_merge   \
     dbsw123      rw_bsw        rw_dat        rw123        \
     zf_bb        zf_cc

# dbsw_tab

bsw_rw: bsw_rw.f90
	$(FC) -o $@ $(link) $< $(DBS) $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

dbound_bsw: dbound_bsw.f90    
	$(FC) -o $@ $(link) $< $(DBS) $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

dbound_tab: dbound_tab.f90
	$(FC) -o $@ $(link) $< $(ZCOM)
	mv $@ $(BINDIR)

dbsr_merge: dbsr_merge.f90        
	$(FC) -o $@ $(link) $< $(DBS) $(ZCONFJJ) $(ZCOM)   -fbacktrace -fcheck=bounds
	mv $@ $(BINDIR)

dbsw123: dbsw123.f90
	$(FC) -o $@ $(link) $< $(DBS) $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

rw123: rw123.f90
	$(FC) -o $@ $(link) $< $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

rw_bsw: rw_bsw.f90
	$(FC) -o $@ $(link) $< $(DBS)  $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

rw_dat: rw_dat.f90    
	$(FC) -o $@ $(link) $< $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

zf_bb: zf_bb.f90
	$(FC) -o $@ $(link) $< $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

zf_cc: zf_cc.f90
	$(FC) -o $@ $(link) $< $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

# -fbacktrace -fcheck=bounds 
