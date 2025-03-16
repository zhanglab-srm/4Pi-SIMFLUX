function [stepnum]=getstepnum_v2(filename)
% inded=find(filename=='.',1,'last');
% indst=find(filename=='_',1,'last');
ind=find(filename=='_');
stepnum=str2double(filename(ind(2)+1:ind(3)-1));