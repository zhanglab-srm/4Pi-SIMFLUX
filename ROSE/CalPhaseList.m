function result = CalPhaseList(intlist, phase0)
%
%result: [phase1 phase2 amp1 amp2 offset1 offset2 md1 md2]
%           1      2     3     4     5      6      7   8
    if nargin<2
        phase0 = 120/180*pi;
    end
    intlen = size(intlist,1);
    result = zeros(intlen,8);
    for m=1:intlen
        [phase1, amp1, offset1] = MyCalPhase3_new(intlist(m,1:3), phase0);
        [phase2, amp2, offset2] = MyCalPhase3_new(intlist(m,4:6), phase0);
        result(m,1) = phase1;
        result(m,3) = amp1;
        result(m,5) = offset1;
        result(m,2) = phase2;
        result(m,4) = amp2;
        result(m,6) = offset2;
    end
    
    result(:,7) = result(:,3)./result(:,5);
    result(:,8) = result(:,4)./result(:,6);
end