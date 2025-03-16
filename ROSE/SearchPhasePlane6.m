function [lp_optimized, result] = SearchPhasePlane6(x,y,p2,samplerate, searchrange)
    %% sample
    if nargin < 5
        searchrange = [-3.2 3.2];
    end
    if nargin <4
        samplerate = 1;
    end
    if samplerate > 1
        samplerate = min(1, samplerate/length(x));
    end
    
    mask = rand(size(x)) <= samplerate;
    xs = x(mask);
    ys = y(mask);
    p2s = p2(mask);
    searchstart = min(searchrange);
    searchend = max(searchrange);
    %% search for parameters
%     lplist1 = searchstart:0.005:searchend;
%     lplist2 = 3.6:0.02:4.0; % 普通照明 3.5:0.02:3.7
% %     lplist2 = 3.9:0.02:4.2; %TIRF照明


   lplist1 = searchstart:0.005:searchend;
    lplist2 = 3.2:0.02:4.5; % 普通照明 3.5:0.02:3.7
%     lplist2 = 3.9:0.02:4.2; %TIRF照明


    [aa, bb] = meshgrid(lplist1, lplist2);
    resnormmap = zeros(size(aa));
    
    lplist = zeros(length(aa(:)), 3);
    %options = optimset('Display','off','MaxFunEvals',100,'MaxIter',5,'TolFun',1e-4,'LargeScale','off');
%     parfor_progress(length(aa(:)));
    parfor m=1:length(aa(:))
        [lp, resnormx] = FitPhasePlane(xs,ys,p2s, [aa(m) bb(m) 0], 1);
        resnormmap(m) = resnormx;
        lplist(m,:) = lp;
%         [lpx,resnormx,residualx,exitflagx]=lsqcurvefit(@(xp, xdata)CalPhaseSin(xp, xs, ys), ...
%         [aa(m) bb(m) 0],xs,[sin(p2s) cos(p2s)],[],[],options);
%         resnormmap(m) = resnormx;
%         lplist(m,:) = lpx;
%         parfor_progress();
    end
%     parfor_progress(0);

    %% get parameters under 500
    mask = find(resnormmap <(min(resnormmap(:))+1));
%     alist = aa(mask);
%     blist = bb(mask);
    lplist = lplist(mask, :);
    options = optimset('Display','off','MaxFunEvals',1000,'MaxIter',100,'TolFun',1e-5,'LargeScale','on');
%     for m=1:length(mask)
    [lpx,resnormx,residualx,exitflagx]=lsqcurvefit(@(xp, xdata)CalPhaseSin(xp, x, y), ...
    [mean(lplist(:,1)) mean(lplist(:,2)) 0],x,[sin(p2) cos(p2)],[],[],options);
%     lplist(m,:) =  lpx;
%     end

    lp_optimized = [lpx(1), lpx(2), checkPhase(lpx(3))];
    result = CalPhase(lp_optimized, x, y);
    %% display
    figure()
    surf(aa,bb,resnormmap);
end