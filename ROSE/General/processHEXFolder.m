function processHEXFolder(dirstring)
    filelist = dir(dirstring);
    t = strfind(dirstring, '\');
    if length(t) >0
        folder = dirstring(1:t(end));
    else
        folder = '';
    end
    filenum = length(filelist);
    for m=1:filenum
        tfile = filelist(m);
        fullpath = [folder tfile.name];
        if tfile.isdir
            continue;
        end
        fprintf('Processing file: %s (%d/%d)\n',fullpath,m,filenum);
        
        processHEXFile(fullpath);
    end
end