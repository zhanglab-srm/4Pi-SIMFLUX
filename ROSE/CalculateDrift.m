function [dx dy distbuf] = CalculateDrift(x1, y1, x2, y2)
maxdist =  4;
plen1 = length(x1);
plen2 = length(x2);
distbuf = zeros(round(plen1*plen2/10), 2);
datalen = 0;

for m=1:plen1
    tx1 = x1(m);
    ty1 = y1(m);
    for n=1:plen2
        tx2 = x2(n);
        ty2 = y2(n);
        if abs(tx1-tx2) < maxdist && abs(ty1-ty2)<maxdist
            datalen=datalen+1;
            distbuf(datalen,1) = tx2-tx1;
            distbuf(datalen,2) = ty2-ty1;
        end
    end
end
distbuf = distbuf(1:datalen,:);
%% Calculate Distance
tx = distbuf(:,1);
ty = distbuf(:,2);
[ax bx] = hist(tx, 40);
[tmean, tstd] = normfit(tx);
[x0x y0x widthx ryx] = GaussianFit1D(bx,ax,tstd);

[ay by] = hist(ty, 40);
[tmean, tstd] = normfit(ty);
[x0y y0y widthy ryy] = GaussianFit1D(by,ay,tstd);

dx = x0x;
dy = x0y;