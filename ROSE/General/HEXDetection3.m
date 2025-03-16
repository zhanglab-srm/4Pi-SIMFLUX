function [result] = HEXDetection3(imgbuf, threshold, detectionsize, disk_radius)
%single molecule detection of image stack
%modified for run in cluster
%[result] = HEXDetection3(imgbuf, threshold, detectionsize, disk_radius)

    if nargin < 4
        disk_radius = 4;
    end
    if nargin < 3
        detectionsize = 15;
    end
    if nargin < 2
        threshold = 4;
    end
    
    maxpacklen = 200;
    minpacklen = 20;
    defpacknum = 200;%per worker
    reqsize = 20;%req buffer size
    p=gcp();
    workernum = p.NumWorkers;
    defpacknum = defpacknum*workernum;
    
    imglen = size(imgbuf,3);
    packlen = min(ceil(imglen/defpacknum), maxpacklen);
    packlen = max(packlen, minpacklen);
    
    packnum = ceil(imglen/packlen);
    result = cell(imglen, 1);
    fprintf('Detection pack Info:\npacklen = %d\npacknum = %d\n', packlen, packnum);
    h = waitbar(0,'Detecting..');
    if packnum <= reqsize   %packlen < reqsize, send request at one time
        fslist = zeros(packnum,1);
        felist = zeros(packnum,1);
        for idx = 1:packnum-1
            fs = (idx-1)*packlen+1;
            fe = idx*packlen;
            f(idx) = parfeval(p,@HEXDetection_st,1,imgbuf(:,:,fs:fe), threshold, detectionsize, disk_radius); 
            fslist(idx) = fs;
            felist(idx) = fe;
        end
        idx = packnum;
        fs = (idx-1)*packlen+1;
        fe = imglen;
        f(idx) = parfeval(p,@HEXDetection_st,1,imgbuf(:,:,fs:fe), threshold, detectionsize, disk_radius); 
        fslist(idx) = fs;
        felist(idx) = fe;
        % Collect the results as they become available.
        for idx = 1:packnum
          % fetchNext blocks until next results are available.
          [completedIdx,tresult] = fetchNext(f);
          fs = fslist(completedIdx);
          fe = felist(completedIdx);
          result(fs:fe,:) = tresult;
          waitbar(idx/packnum, h, sprintf('Detecting.. %02.1f%%(%d/%d)',idx/packnum*100,idx, packnum));
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
            f(idx) = parfeval(p,@HEXDetection_st,1,imgbuf(:,:,fs:fe), threshold, detectionsize, disk_radius); 
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
            f(completedIdx) = parfeval(p,@HEXDetection_st,1,imgbuf(:,:,fs:fe), threshold, detectionsize, disk_radius); 
            fslist(completedIdx) = fs;
            felist(completedIdx) = fe;
            sendcnt = sendcnt + 1;
%             disp(num2str(fetchcnt));
            waitbar(fetchcnt/packnum, h, sprintf('Detecting.. %02.1f%%(%d/%d)',fetchcnt/packnum*100,fetchcnt, packnum));
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
            waitbar(fetchcnt/packnum, h, sprintf('Detecting.. %02.1f%%(%d/%d)',fetchcnt/packnum*100,fetchcnt, packnum));
        end
    end
    
    close(h);
end