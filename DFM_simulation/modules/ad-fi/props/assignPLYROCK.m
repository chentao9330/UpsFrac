function f = assignPLYROCK(f, plyrock, reg)
ntsfun = numel(reg.SATINX);
if ntsfun == 1
    pvtnum = 1;
else
    pvtnum = reg.PVTNUM;
end
f.dps    = plyrock(pvtnum, 1);
f.rrf    = plyrock(pvtnum, 2);
f.rhoR   = plyrock(pvtnum, 3);
f.adsInx = plyrock(pvtnum, 4);
f.adsMax   = plyrock(pvtnum, 5);
end
