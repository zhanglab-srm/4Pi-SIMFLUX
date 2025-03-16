function result = Calstd(cnt,bin)
    cnt=cnt-min(cnt);
    meany = sum(bin.*cnt)./sum(cnt);
    meany2 = sum(bin.^2.*cnt)./sum(cnt);
    
    result = sqrt(meany2 - meany^2);
end