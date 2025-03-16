function [detectResults] = DWTParticalDetection(img, threshold, dettype, windowWidth)
%     kernel{1}=[1/16,1/4,3/8,1/4,1/16];
%     kernel{2}=[1/16,0,1/4,0,3/8,0,1/4,0,1/16];
%     kernel{3}=[1/16,0,0,0,1/4,0,0,0,3/8,0,0,0,1/4,0,0,0,1/16];
%     [tW2 tW3] = det_DWT_cpu_pad_float(single(img),single(kernel{1}), single(kernel{2}),single(kernel{3}));
    kernel_compressed = [3/8 1/4 1/16];
%     [tW2 tW3] = det_DWT_newconv2(single(img),single(kernel_compressed));

%     W2 = det_Thresh(tW2, threshold);
%     W3 = det_Thresh(tW3, threshold);
%     W2 = det_DWT_thresh(tW2, threshold);
%     W3 = det_DWT_thresh(tW3, threshold);
    
    samplewidth = round((windowWidth -1)/2);
    switch(dettype)
        case 1
            Wresult = det_DWT_newconv2_threshold(single(img), single(kernel_compressed), single(threshold), 'W2');
            detectResults = det_DWT_Findparticles_float(Wresult, samplewidth);
        case 2
            Wresult = det_DWT_newconv2_threshold(single(img), single(kernel_compressed), single(threshold), 'W3');
            detectResults = det_DWT_Findparticles_float(Wresult, samplewidth);
        case 3
            Wresult = det_DWT_newconv2_threshold(single(img), single(kernel_compressed), single(threshold), 'W2+W3');
            detectResults = det_DWT_Findparticles_float(Wresult, samplewidth);
        case 4
            Wresult = det_DWT_newconv2_threshold(single(img), single(kernel_compressed), single(threshold), 'W2*W3');
            detectResults = det_DWT_Findparticles_float(Wresult, samplewidth);
    end
end