function result = CalPhaseListLite(intlist, phase0)
%
%result: [phase1  amp1  offset1  md1 ]
%           1      2       3      4  
    if nargin<2
        phase0 = 120/180*pi;
    end
    intlen = size(intlist,1);
    result = zeros(intlen,4);
    for m=1:intlen
        [phase1, amp1, offset1] = MyCalPhase3_new(intlist(m,1:3), phase0);
%         [phase2, amp2, offset2] = MyCalPhase3_new(intlist(m,4:6), phase0);
        result(m,1) = phase1;
        result(m,2) = amp1;
        result(m,3) = offset1;
%         result(m,2) = phase2;
%         result(m,4) = amp2;
%         result(m,6) = offset2;
    end
    
    result(:,4) = result(:,2)./result(:,3);
%     result(:,8) = result(:,4)./result(:,6);
end