function save_hdf5(fileName, data, datalen)%data as structure
%save_hdf5(fileName, data, datalen)
%save data into hdf5 file
%data:
%x, y, frame, int/photon, sx, sy
%

%% hdf5 file test
% data in structure
tdata = [];
datalen = len(data.x);
tdata.frame = zeros(1,datalen,'uint32');
tdata.x=zeros(1,datalen,'single');
tdata.y=zeros(1,datalen,'single');
tdata.photons=zeros(1,datalen,'single');
tdata.sx=ones(1,datalen,'single').*0.05;
tdata.sy=ones(1,datalen,'single').*0.05;
tdata.bg=ones(1,datalen,'single');
tdata.lpx=ones(1,datalen,'single').*0.05;
tdata.lpy=ones(1,datalen,'single').*0.05;
tdata.net_gradient=ones(1,datalen,'single');
tdata.likelihood=ones(1,datalen,'single');
tdata.iterations=ones(1,datalen,'uint32');

if isfield(data,'x')
    tdata.x = single(data.x);
%     datalen = len(data.x);
end

if isfield(data,'y')
    tdata.y = single(data.y);
%     datalen = len(data.y);
end

if isfield(data,'frame')
    tdata.frame = uint32(data.frame);
%     datalen = len(data.frame);
end

if isfield(data,'photon')
    tdata.photons = single(data.photon);
%     datalen = len(data.photon);
end

if isfield(data,'int')
    tdata.photons = single(data.int);
%     datalen = len(data.int);
end

if isfield(data,'sx')
    tdata.sx = single(data.sx);
%     datalen = len(data.sx);
end
if isfield(data,'sy')
    tdata.sy = single(data.sy);
%     datalen = len(data.sy);
end

% fileName       = 'test.hdf5';
DATASET        = 'locs';
DIM0           = datalen;

dims = DIM0;

%% create file
%
file = H5F.create (fileName, 'H5F_ACC_TRUNC',...
    'H5P_DEFAULT', 'H5P_DEFAULT');

%
%Create the required data types
%
uintType   =H5T.copy('H5T_NATIVE_UINT32');
sz(1)     =H5T.get_size(uintType);
singleType=H5T.copy('H5T_NATIVE_FLOAT');
sz(2:11)     =H5T.get_size(singleType);
intType   =H5T.copy('H5T_NATIVE_INT32');
sz(12)     =H5T.get_size(intType);


%
% Computer the offsets to each field. The first offset is always zero.
%
offset(1)=0;
offset(2:12)=cumsum(sz(1:11));

%
% Create the compound datatype for memory.
%
memtype = H5T.create ('H5T_COMPOUND', sum(sz));
H5T.insert (memtype,...
    'frame',offset(1),uintType);
H5T.insert (memtype,...
    'x',offset(2), singleType);
H5T.insert (memtype,...
    'y',offset(3), singleType);
H5T.insert (memtype,...
    'photons',offset(4), singleType);
H5T.insert (memtype,...
    'sx',offset(5), singleType);
H5T.insert (memtype,...
    'sy',offset(6), singleType);
H5T.insert (memtype,...
    'bg',offset(7), singleType);
H5T.insert (memtype,...
    'lpx',offset(8), singleType);
H5T.insert (memtype,...
    'lpy',offset(9), singleType);
H5T.insert (memtype,...
    'net_gradient',offset(10), singleType);
H5T.insert (memtype,...
    'likelihood',offset(11), singleType);
H5T.insert (memtype,...
    'iterations',offset(12), intType);

%
% Create the compound datatype for the file.  Because the standard
% types we are using for the file may have different sizes than
% the corresponding native types, we must manually calculate the
% offset of each member.
%
filetype = H5T.create ('H5T_COMPOUND', sum(sz));
H5T.insert (filetype,...
    'frame',offset(1),uintType);
H5T.insert (filetype,...
    'x',offset(2), singleType);
H5T.insert (filetype,...
    'y',offset(3), singleType);
H5T.insert (filetype,...
    'photons',offset(4), singleType);
H5T.insert (filetype,...
    'sx',offset(5), singleType);
H5T.insert (filetype,...
    'sy',offset(6), singleType);
H5T.insert (filetype,...
    'bg',offset(7), singleType);
H5T.insert (filetype,...
    'lpx',offset(8), singleType);
H5T.insert (filetype,...
    'lpy',offset(9), singleType);
H5T.insert (filetype,...
    'net_gradient',offset(10), singleType);
H5T.insert (filetype,...
    'likelihood',offset(11), singleType);
H5T.insert (filetype,...
    'iterations',offset(12), intType);


%
% Create dataspace.  Setting maximum size to [] sets the maximum
% size to be the current size.
%
space = H5S.create_simple (1,fliplr( dims), []);

%
% Create the dataset and write the compound data to it.
%
dset = H5D.create (file, DATASET, filetype, space, 'H5P_DEFAULT');
H5D.write (dset, memtype, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', tdata);

%
% Close and release resources.
%
H5D.close (dset);
H5S.close (space);
H5T.close (filetype);
H5F.close (file);


%% backup
% type_id = H5T.create('H5T_COMPOUND',12);
% H5T.insert(type_id,'frame',0,'H5T_NATIVE_UINT32');
% H5T.insert(type_id,'x',4,'H5T_NATIVE_FLOAT');
% H5T.insert(type_id,'y',8,'H5T_NATIVE_FLOAT');
% 
% h5write('e:\myfile.h5', '/locs', tdata);