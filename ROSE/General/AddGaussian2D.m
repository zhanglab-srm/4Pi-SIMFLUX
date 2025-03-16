function img = AddGaussian2D(img, x, y, height, std)
    s = size(img);
    [xx yy] = meshgrid(1:s(2), 1:s(1));
    
    result = fun_Gaussian(xx,yy,x,y,std);
    img = img + result.*height;

end

function result = fun_Gaussian(x,y,x0,y0,std)
    result = exp(-((x-x0).^2 + (y-y0).^2)/(2*std^2))/(2*pi*std^2);
end