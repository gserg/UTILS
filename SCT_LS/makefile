#          utils SCT_LS

include ${HOME}/DEFS_03

all: add_farm            add_stgf          add_stgfb         adf04          \
     beb_ion             delete_data       eff_cstrn         get1tdonejj    \
     gettdonejk          gettdoneLS        kma_om            kma_phase      \
     kma_phase_lsp       oma_omb           sec_cascb         sec_dif_JK     \
     sec_dif_JK_ampl     sec_dif_LS        sec_dif_LS_tdone  sec_omn        \
     sec_par             sec_parb          sec_sums          sec_top        \
     sec_top_CBE         sec_totalb        tma_tmb           tmb_omb 
     

add_farm: add_farm.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

add_stgf: add_stgf.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)                                                                       

add_stgfb: add_stgfb.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

adf04: adf04.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

beb_ion: beb_ion.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

delete_data: delete_data.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

eff_cstrn: eff_cstrn.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

get1tdonejj: get1tdonejj.f90                               
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

gettdonejk: gettdonejk.f90                               
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

gettdoneLS: gettdoneLS.f90                                
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

kma_om: kma_om.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

kma_phase: kma_phase.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)  $(LAPACK)
	mv $@ $(BINDIR)

kma_phase_lsp: kma_phase_lsp.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)  $(LAPACK)
	mv $@ $(BINDIR)

oma_omb: oma_omb.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

sec_cascb: sec_cascb.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

sec_dif_JK: sec_dif_JK.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

sec_dif_JK_ampl: sec_dif_JK_ampl.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

sec_dif_LS: sec_dif_LS.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

sec_dif_LS_tdone: sec_dif_LS_tdone.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

sec_omn: sec_omn.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

sec_par: sec_par.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

sec_parb: sec_parb.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

sec_sums: sec_sums.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

sec_top: sec_top.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

sec_top_CBE: sec_top_CBE.f90
	$(FC) -o $@ $(link) $< fdip_org.f $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

sec_totalb: sec_totalb.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

tma_tmb: tma_tmb.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

tmb_omb: tmb_omb.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

#  -fbacktrace  -fcheck=bounds  

