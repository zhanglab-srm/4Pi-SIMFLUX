function [dx, dy, distbuf] = CalculateDrift_CPU(x1, y1, x2, y2)
maxdist =  4;
binsize = 101;
% plen1 = length(x1);
% plen2 = length(x2);
% distbuf = zeros(round(plen1*plen2/100), 2);
% datalen = 0;
% 
% for m=1:plen1
%     tx1 = x1(m);
%     ty1 = y1(m);
%     for n=1:plen2
%         tx2 = x2(n);
%         ty2 = y2(n);
%         if abs(tx1-tx2) < maxdist && abs(ty1-ty2)<maxdist
%             datalen=datalen+1;
%             distbuf(datalen,1) = tx2-tx1;
%             distbuf(datalen,2) = ty2-ty1;
%         end
%     end
% end
% distbuf = distbuf(1:datalen,:);
[histx, histy, distbin] = CalPointsetDist_mex(x1,y1,x2,y2,maxdist,binsize);
%% Calculate Distance
tstd = Calstd(histx, distbin);
[x0x y0x widthx ryx] = GaussianFit1D(distbin,histx,tstd);


tstd = Calstd(histy, distbin);
[x0y y0y widthy ryy] = GaussianFit1D(distbin,histy,tstd);

dx = x0x;
dy = x0y;
distbuf = [distbin, histx, histy];
end

