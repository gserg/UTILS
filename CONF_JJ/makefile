#       utilities DBSR

include ${HOME}/DEFS_03

all: cm_j         coef_jj    genjconf    genjterm     \
     jcfile       jfile      jj_list     mult_jj      order_jj     

cm_j: cm_j.f90
	$(FC) -o $@ $(link)  $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

coef_jj: coef_jj.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

genjconf: genjconf.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

genjterm: genjterm.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

jcfile: jcfile.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

jfile: jfile.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

jj_list: jj_list.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

mult_jj: mult_jj.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)

order_jj: order_jj.f90
	$(FC) -o $@ $(link) $<  $(ZCONFJJ) $(ZCOM) 
	mv $@ $(BINDIR)


# -fbacktrace -fcheck=bounds 
