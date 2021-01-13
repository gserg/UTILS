#          utils CONF_LS

include ${HOME}/DEFS_03

all:  cfile     cfile_all    ci_list     del_c       del_core  \
      lfile     merge        order_c     order_cc              \
      zf_tab    zgenconf     zgenterm    

cfile: cfile.f90
	$(FC) -o $@ $(link) $<  
	mv $@ $(BINDIR)

cfile_all: cfile_all.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

ci_list: ci_list.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

del_c: del_c.f90
	$(FC) -o $@ $(link) $<  
	mv $@ $(BINDIR)

del_core: del_core.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

lfile: lfile.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

merge: merge.f90
	$(FC) -o $@ $(link) $<  $(RADIAL) $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

mfile: mfile.f90
	$(FC) -o $@ $(link) $<  $(RADIAL) $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

order_c: order_c.f90
	$(FC) -o $@ $(link) $<  
	mv $@ $(BINDIR)

order_cc: order_cc.f90
	$(FC) -o $@ $(link) $<  
	mv $@ $(BINDIR)

zf_tab: zf_tab.f90
	$(FC) -o $@ $(link) $<  $(ZCOM)
	mv $@ $(BINDIR)

zf_tab_E2: zf_tab_E2.f90
	$(FC) -o $@ $(link) $<  $(ZCOM)
	mv $@ $(BINDIR)

zf_tab_exp: zf_tab_exp.f90
	$(FC) -o $@ $(link) $<  $(ZCOM)
	mv $@ $(BINDIR)

zf_tab_jj: zf_tab_jj.f90
	$(FC) -o $@ $(link) $<  $(ZCOM)
	mv $@ $(BINDIR)

zgenconf: zgenconf.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)

zgenterm: zgenterm.f90
	$(FC) -o $@ $(link) $<  $(ZCONFLS) $(ZCOM)
	mv $@ $(BINDIR)
