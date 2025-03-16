function poslist = HEXDetection_single(img, threshold, detectionsize, disk_radius)
    if nargin < 4
        disk_radius = 4;
    end
    if nargin < 3
        detectionsize = 15;
    end
    if nargin < 2
        threshold = 4;
    end
    % Filter with Gaussian blur
    
    h = GenHEX(11,[0,0,2.0,1,0,0,disk_radius],0);
    h = h - mean(h(:));
    
    timg = double(img);
    t1 = medfilt2(timg, [detectionsize detectionsize]);
    t2 = timg - t1;
    t3 = filter2(h, t2);
    mask1 = (t3>(mean(t3(:)) + std(t3(:))*threshold));
    fImg=locmax2d(t3, [detectionsize detectionsize]);

    t4 = fImg>0 & mask1;

    t = find(t4);
    [ty, tx] = ind2sub(size(t4), t);
    
    poslist = [tx ty];
end