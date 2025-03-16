%% fit six images
%data in timgbuf in saze n*n*6*imglen
%parameter: [x y wx wy int1 int2 int3 int4 int5 int6 bkg1 bkg2 bkg3 bkg4 bkg5 bkg6 exitflag]
%           [1 2  3  4  5    6    7    8    9    10   11    12   13   14   15   16    17]
function result = fit6Gauss4_st(timgbuf_)
    [result, info1, info2] = GaussFit6_CPU(timgbuf_);
end