function [result, fitdata, phasedata] = processBeadsFile(filename, forcefit)
%[result, fitdata, phasedata] = processHeadsFile(filename)
%result: [x y p1 p2 std int md1 md2]
%         1 2 3   4  5   6   7  8  
%fitdata: [x0 y0 std int1 int2 int3 int4 int5 int6 bkg rx ry phase_offset]
%          1  2  3   4    5    6    7     8    9    10 11 12     13
%phasedata: [phase1, amp1, offset1, phase2, amp2, offset2]
%              1      2       3       4      5      6 
if nargin <2
    forcefit = 0;
end
if forcefit ==0 && exist([filename(1:end-4) '_fitdata.mat'], 'file')
    result=[];
    fitdata=[];
    phasedata=[];
    load([filename(1:end-4) '_fitdata.mat'])
    if ~isempty(result) && ~isempty(fitdata) && ~isempty(phasedata)
        return
    end
end
imgbuf = tiffread(filename);
[result, fitdata, phasedata] = fitHEX(imgbuf);
save([filename(1:end-4) '_fitdata.mat'], 'result', 'fitdata', 'phasedata');
end