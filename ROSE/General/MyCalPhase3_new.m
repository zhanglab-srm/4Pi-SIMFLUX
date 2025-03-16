function [phase, amp, offset] = MyCalPhase3_new(intlist, phase0)
    int1 = intlist(1);
    int2 = intlist(2);
    int3 = intlist(3);
    cp = cos(phase0);
    sp = sin(phase0);
    
    offset = (int2*cp - (int1+int3)/2)/(cp-1);
    int1 = int1 - offset;
    int2 = int2 - offset;
    int3 = int3 - offset;
    
    amp = sqrt(int2^2 + (int1-int2*cp)^2/(sp^2));
    
    sx = int2/amp;
    cx = (int1 - int2*cp)/sp/amp;
    
    tp = acos(cx);
    if sx <0
        phase = -real(tp);
    else
        phase = real(tp);
    end
end