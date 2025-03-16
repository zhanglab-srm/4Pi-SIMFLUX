function [phaseBuf, md] = MyCalPhase(intlist, phase0)
    if nargin < 2
        phase0 = 120/180*pi;
    end
    imglen = size(intlist,1);
    phaseBuf = zeros(imglen,6);%phase1, amp1, offset1, phase2, amp2, offset2
    for m=1:imglen
        tintlist = intlist(m, :);
        offset = min(tintlist);
        tintlist = tintlist - offset;
        [phaseBuf(m,1), phaseBuf(m,2), phaseBuf(m,3)] = MyCalPhase3_new(tintlist(1:3), phase0);
        phaseBuf(m,3) = phaseBuf(m,3) + offset;
    end
    md = phaseBuf(:,2)./phaseBuf(:,3);
end