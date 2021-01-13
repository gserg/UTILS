zgenterm cv1.conf 0 0 5 5 4 4 cv1.c
copy cv1.c cfg.inp
nonh
copy 3p6_3d.w wfn.inp
mchf < m.inp
copy summry cv1.s
copy wfn.out cv1.w
ci < ci.inp




