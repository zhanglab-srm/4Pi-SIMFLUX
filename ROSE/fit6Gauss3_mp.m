function [result] = fit6Gauss3_mp(imgbuf)
% fit and calculate phase of image stack
% modified for run in cluster
% [result] = fitHEX(imgbuf)
% data in timgbuf in saze n*n*6*imglen
% result:    [x y wx wy int1 int2 int3 int4 int5 int6 bkg1 bkg2 bkg3 bkg4 bkg5 bkg6 exitflag]
%            [1 2  3  4  5    6    7    8    9    10   11    12   13   14   15   16    17]
    
    maxpacklen = 2000;
    minpacklen = 100;
    defpacknum = 100;%per worker
    reqsize = 20;
    p=gcp();
    workernum = p.NumWorkers;
    defpacknum = defpacknum*workernum;
    
    imglen = size(imgbuf,4);
    packlen = min(ceil(imglen/defpacknum), maxpacklen);
    packlen = max(packlen, minpacklen);
    
    packnum = ceil(imglen/packlen);
    result = zeros(imglen, 17);
    fprintf('packInfo:\npacklen = %d\npacknum = %d\n', packlen, packnum);
    h = waitbar(0,'Fitting..');
    if packnum <= reqsize   %packlen < reqsize, send request at one time
        fslist = zeros(packnum,1);
        felist = zeros(packnum,1);
        for idx = 1:packnum-1
            fs = (idx-1)*packlen+1;
            fe = idx*packlen;
            f(idx) = parfeval(p,@fit6Gauss3_st,1,imgbuf(:,:,:,fs:fe)); 
            fslist(idx) = fs;
            felist(idx) = fe;
        end
        idx = packnum;
        fs = (idx-1)*packlen+1;
        fe = imglen;
        f(idx) = parfeval(p,@fit6Gauss3_st,1,imgbuf(:,:,:,fs:fe));
        fslist(idx) = fs;
        felist(idx) = fe;
        % Collect the results as they become available.
        for idx = 1:packnum
          % fetchNext blocks until next results are available.
          [completedIdx,tresult] = fetchNext(f);
          fs = fslist(completedIdx);
          fe = felist(completedIdx);
          result(fs:fe,:) = tresult;
          waitbar(idx/packnum, h, sprintf('Fitting.. %02.1f%%(%d/%d)',idx/packnum*100,idx, packnum));
        end
    else  %packnum > reqsize, send request WITH fetch
        sendcnt = 0;
        fetchcnt = 0;
        %Step1: fill request buffer
        fslist = zeros(reqsize,1);
        felist = zeros(reqsize,1);
        for idx = 1:reqsize
            fs = (sendcnt)*packlen+1;
            fe = (sendcnt+1)*packlen;
            f(idx) = parfeval(p,@fit6Gauss3_st,1,imgbuf(:,:,:,fs:fe));
            fslist(idx) = fs;
            felist(idx) = fe;
            sendcnt = sendcnt +1;
        end
        %Step2: fetch result and refill request untill all request are send
        while sendcnt < packnum
            % fetchNext blocks until next results are available.
            [completedIdx,tresult] = fetchNext(f);
            fs = fslist(completedIdx);
            fe = felist(completedIdx);
            result(fs:fe,:) = tresult;
            fetchcnt = fetchcnt + 1;
            %send new request
            fs = (sendcnt)*packlen+1;
            fe = min(imglen, (sendcnt+1)*packlen);
            f(completedIdx) = parfeval(p,@fit6Gauss3_st,1,imgbuf(:,:,:,fs:fe)); 
            fslist(completedIdx) = fs;
            felist(completedIdx) = fe;
            sendcnt = sendcnt + 1;
%             disp(num2str(fetchcnt));
            waitbar(fetchcnt/packnum, h, sprintf('Fitting.. %02.1f%%(%d/%d)',fetchcnt/packnum*100,fetchcnt, packnum));
        end
        %Step3: fetch result
        while fetchcnt < packnum
            % fetchNext blocks until next results are available.
            [completedIdx,tresult] = fetchNext(f);
            fs = fslist(completedIdx);
            fe = felist(completedIdx);
            result(fs:fe,:) = tresult;
            fetchcnt = fetchcnt + 1;
%             disp(num2str(fetchcnt));
            waitbar(fetchcnt/packnum, h, sprintf('Fitting.. %02.1f%%(%d/%d)',fetchcnt/packnum*100,fetchcnt, packnum));
        end
    end
    
    close(h);
end