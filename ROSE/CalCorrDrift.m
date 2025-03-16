function [driftx, drifty, dx, dy, xcorr, ycorr] = CalCorrDrift(x,y,frame,minpts,minframes)
%Calculate drift data with correlation methods
%[dx, dy, xcorr, ycorr] = CalCorrDrift(x,y,frame,minpts,minframes)
%
%
%
    %% data in x, y, frame
    Xr = x;
    Yr = y;
    % frame = fposlist(:,5);
    if nargin <4
        minpts = 5000;
    end
    if nargin <5
        minframes = 50;
    end

    plen = length(Xr);
    framelen = max(frame) - min(frame)+1;
    pnum = minpts;
    if framelen/(plen/pnum) <minframes  %increase pnum
        pnum = round(plen*minframes/framelen/2)*2;
    end
    
    step = pnum/2;
    stepcnt = ceil(plen/step);

    driftbuf = zeros(ceil(plen./step), 4);
    driftlen = 0;
    start = 1;
    x1 = Xr(start:(start+pnum));
    y1 = Yr(start:(start+pnum));
    f1 = frame(start:(start+pnum));
    rawdata = cell(ceil(plen./step), 1);

    h=waitbar(0,'Drift processing..');
    while(start+pnum+step < plen)

        x2 = Xr((start+step):(start+pnum+step));
        y2 = Yr((start+step):(start+pnum+step));
        f2 = frame((start+step):(start+pnum+step));
        [dx dy distbuf] = CalculateDrift(x1, y1, x2, y2);
        driftlen = driftlen +1;
        driftbuf(driftlen,1) = dx;
        driftbuf(driftlen,2) = dy;
        driftbuf(driftlen,3) = mean(f1);
        driftbuf(driftlen,4) = mean(f2);
        start = start + step;
        rawdata{driftlen} = distbuf;
        waitbar(driftlen/stepcnt,h);
    end
    close(h);
    disp(sprintf('Correlation Drift Step: %d',driftlen));

    driftbuf = driftbuf(1:driftlen, :);
    rawdata = rawdata(1:driftlen);

    %% drift data
    framelist = 1:max(frame);
    driftx = interp1(driftbuf(:,4), driftbuf(:,1),framelist,'linear','extrap');
    drifty = interp1(driftbuf(:,4), driftbuf(:,2),framelist,'linear','extrap');

    % [lp, driftx_fit, exitflag] = FitLine(framelist, driftx);
    % [lp, drifty_fit, exitflag] = FitLine(framelist, drifty);

    dx = driftx(frame);
    dy = drifty(frame);
    dx=dx';
    dy=dy';
    %% drift correction
    xcorr = x - dx;
    ycorr = y - dy;
end