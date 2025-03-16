function img = PointSet2Image(xlist, ylist, zoom, rect)
%img = PointSet2Image(xlist, ylist, zoom, rect)
    if nargin <4
        minx = floor(min(xlist));
        maxx = ceil(max(xlist));
        miny = floor(min(ylist));
        maxy = ceil(max(ylist));
    else
        minx = rect(1);
        maxx = rect(2);
        miny = rect(3);
        maxy = rect(4);
    end
    
    mask = xlist>=minx & xlist<=maxx & ylist>=miny & ylist<=maxy;
    xlist= xlist(mask);
    ylist= ylist(mask);
    
%     img = zeros((maxy-miny+2)*zoom, (maxx-minx+2)*zoom, 'int16');

    img = zeros((maxy-miny+1)*zoom, (maxx-minx+1)*zoom, 'int16');
    
    xlist =round((xlist - minx + 1)*zoom);
    ylist =round((ylist - miny + 1)*zoom);
    for m=1:length(xlist)
        img(ylist(m), xlist(m)) = img(ylist(m), xlist(m)) +1;
    end
end