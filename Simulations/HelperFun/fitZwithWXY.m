function retz = fitZwithWXY(wxlist, wylist, wxparam, wyparam, zrange, zstep)
if nargin<5
    zrange=[-200 1200];
end
if nargin<6
    zstep=2;
end
zlist = zrange(1):zstep:zrange(2);
twxlist = polyval(wxparam, zlist);
twylist = polyval(wyparam, zlist);
retz = fitz_wxy_mex(zlist, twxlist, twylist, wxlist, wylist);