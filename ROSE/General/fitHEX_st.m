function [result, fitdata, phasedata] = fitHEX_st(imgbuf, r)
%Single Threaded version of fitHEX
%fit and calculate phase of image stack
%[result fitdata phasedata] = fitHEX(imgbuf)
%result: [x y p1 p2 std int md1 md2]
%         1 2 3   4  5   6   7   8   
%fitdata: [x0 y0 std int1 int2 int3 int4 int5 int6 bkg rx ry phase_offset]
%          1  2  3   4    5    6    7     8    9    10 11 12     13
%phasedata: [phase1, amp1, offset1, phase2, amp2, offset2]
%              1      2       3       4      5      6 
    if nargin <2
        r=4;
    end
    imgbuf = double(imgbuf);
%     s = size(imgbuf);
    imglen = size(imgbuf,3);
    resultbuf = zeros(imglen, 13);
    for m=1:imglen
        timg = imgbuf(:,:,m);
        ret = fitGaussian_Hex_Oval(timg, r);
        resultbuf(m,:) = ret;
    end

    %% Phase
    phase0 = 120/180*pi;
    phaseBuf = zeros(imglen,6);%phase1, amp1, offset1, phase2, amp2, offset2
    for m=1:imglen
        intlist = resultbuf(m, 4:9);
        [phaseBuf(m,1), phaseBuf(m,2), phaseBuf(m,3)] = MyCalPhase3_new(intlist(1:3), phase0);
        [phaseBuf(m,4), phaseBuf(m,5), phaseBuf(m,6)] = MyCalPhase3_new(intlist(4:6), phase0);
    end
    md1 = phaseBuf(:,2)./phaseBuf(:,3);
    md2 = phaseBuf(:,5)./phaseBuf(:,6);
    
    %% make results
    %result: [x y p1 p2 std int md1 md2]
    %         1 2 3   4  5   6  7    8  
    %fitdata: [x0 y0 std int1 int2 int3 int4 int5 int6 bkg rx ry phase_offset]
    %          1  2  3   4    5    6    7     8    9    10 11 12     13
    %phasedata: [phase1, amp1, offset1, phase2, amp2, offset2]
    %              1      2       3       4      5      6 
    fitdata = resultbuf;
    phasedata = phaseBuf;
    
    result = zeros(imglen, 8);
    result(:,1) = fitdata(:,1);
    result(:,2) = fitdata(:,2);
    result(:,3) = phasedata(:,1);
    result(:,4) = phasedata(:,4);
    result(:,5) = fitdata(:,3);
    result(:,6) = sum(fitdata(:,4:9),2);
    result(:,7) = md1;
    result(:,8) = md2;
end