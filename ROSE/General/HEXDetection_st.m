function result = HEXDetection_st(imgbuf, threshold, detectionsize, disk_radius)   
    imglen = size(imgbuf,3);
    result = cell(imglen,1);
    for m=1:imglen
        result{m} = HEXDetection_single(imgbuf(:,:,m), threshold, detectionsize, disk_radius);
    end
end