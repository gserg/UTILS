#          utils SCT_JJ

include ${HOME}/DEFS_03

all:  f_values_jj     get1tdonejj       gettdonejj   \
      sec_dif_JJ      sec_dif_JJ_ampl   sum_hh_jj    \
      targ_jj         targ_jj_exp

f_values_jj: f_values_jj.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

get1tdonejj: get1tdonejj.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)
	rm *.mod

gettdonejj: gettdonejj.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)
	rm *.mod

sec_dif_JJ: sec_dif_JJ.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

sec_dif_JJ_ampl: sec_dif_JJ_ampl.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

sum_hh_jj: sum_hh_jj.f90
	$(FC) -o $@ $(link) $<  $(ZCOM)
	mv $@ $(BINDIR)

targ_jj: targ_jj.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

targ_jj_exp: targ_jj_exp.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM)
	mv $@ $(BINDIR)

# -fbacktrace -fcheck=bounds