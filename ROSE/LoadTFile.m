function [T, poslist_source, poslist_target] = LoadTFile(filepath)
%% load tranform file
% filepath = 'E:\ROSE\20181012_2CCD_Registration\2.txt';
fp = fopen(filepath);
buf = fread(fp);
fclose(fp);

buf = char(buf');

plist = strfind(buf,sprintf('\n'));
p1 = strfind(buf,'Refined source landmarks');
p2 = strfind(buf,'Target landmarks');

plist1 = plist(plist>p1);
plist2 = plist(plist>p2);

str1 = buf(plist1(1)+1 : plist1(2)-1);
str2 = buf(plist1(2)+1 : plist1(3)-1);
str3 = buf(plist1(3)+1 : plist1(4)-1);
str4 = buf(plist2(1)+1 : plist2(2)-1);
str5 = buf(plist2(2)+1 : plist2(3)-1);
str6 = buf(plist2(3)+1 : plist2(4)-1);

% convert to double
pos1 = (str2double(strsplit(str1))) +1;
pos2 = (str2double(strsplit(str2))) +1;
pos3 = (str2double(strsplit(str3))) +1;
pos4 = (str2double(strsplit(str4))) +1;
pos5 = (str2double(strsplit(str5))) +1;
pos6 = (str2double(strsplit(str6))) +1;

poslist_source = cat(1, pos1,pos2,pos3);
poslist_target = cat(1, pos4,pos5,pos6);

%% calculate the tranform matrix T
A = cat(2, poslist_source, ones(size(poslist_source,1),1));
B = cat(2, poslist_target, ones(size(poslist_target,1),1));

T = (A'*A)\(A'*B);