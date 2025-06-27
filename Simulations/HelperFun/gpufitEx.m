function [parameters, states, chi_squares, n_iterations, time]...
    = gpufitEx(data, model_id, initial_parameters, varargin)
% [parameters, states, chi_squares, n_iterations, time]...
%     = gpufitEx(data, model_id, initial_parameters, varargin)
% Wrapper around the Gpufit mex file.
%
%
% Default values as specified
%
% This Ex function contains:
% * Auto data shape conversion
% * Auto data type conversion
% * Auto initial parameter expansion
% * Returned parameters is converted to double
%

% Required parameters:
% data, model_id, initial_parameters
% 
% Optional parameters compared with parameters in gpufit:
% weight    : weights
% tol       : tolerance
% max_iter  : max_n_iterations
% par_mask  : parameters_to_fit
% est_id    : estimator_id
% user_info : user_info
%
%% parameter assign

% number of input parameter (variable)
if nargin < 3
    error('GPUfit: Not enough parameters');
end
weights = [];
tolerance = [];
max_n_iterations = [];
parameters_to_fit = [];
estimator_id = [];
user_info = [];

if mod(len(varargin),2) >0
    error('GPUfit: optional parameter length error');
end
% assign optional parameters
for m=1:2:len(varargin)
    cparname = varargin{m};
    switch(cparname)
        case 'weight' 
            weights = varargin{m+1};
        case 'tol'
            tolerance = varargin{m+1};
        case 'max_iter'
            max_n_iterations = varargin{m+1};
        case 'par_mask'
            parameters_to_fit = varargin{m+1};
        case 'est_id'
            estimator_id = varargin{m+1};
        case 'user_info'
            user_info = varargin{m+1};
    end
end

%% size checks
% data is 2D and read number of points and fits
data_size = size(data);  % data: [7,7,6,3734]  x范围，y范围，6张图，3734个单分子点

if length(data_size) >2
    %data is not two-dimensional, converted to 2 demensional
    disp('gpufitEx: Notice: data is not 2D in GPUfit, converting to 2-D shape');
    n_points = 1;
    for m=1:(length(data_size)-1)
        n_points = n_points*data_size(m);
    end              
    n_fits = data_size(end);
    data = reshape(data, [n_points, n_fits]);
    data_size = size(data);
end
n_points = data_size(1);
n_fits = data_size(2);

% consistency with weights (if given)
if ~isempty(weights)
    assert(isequal(data_size, size(weights)), 'Dimension mismatch between data and weights')
end

% initial parameters is 2D and read number of parameters
initial_parameters_size = size(initial_parameters);
assert(length(initial_parameters_size) == 2, 'initial_parameters is not two-dimensional');
n_parameters = initial_parameters_size(1);

%expand initial_parameters if initial_parameters contains only one column
if initial_parameters_size(2)==1
    initial_parameters = repmat(initial_parameters, 1, n_fits);
    initial_parameters_size = size(initial_parameters);
end
assert(n_fits == initial_parameters_size(2), 'Dimension mismatch in number of fits between data and initial_parameters');

% consistency with parameters_to_fit (if given)
if ~isempty(parameters_to_fit)
    assert(size(parameters_to_fit, 1) == n_parameters, 'Dimension mismatch in number of parameters between initial_parameters and parameters_to_fit');
end

%% default values

% tolerance
if isempty(tolerance)
    tolerance = 1e-6;
end

% max_n_iterations
if isempty(max_n_iterations)
    max_n_iterations = 50;
end

% estimator_id
if isempty(estimator_id)
    estimator_id = EstimatorID.LSE;
end

% parameters_to_fit
if isempty(parameters_to_fit)
    parameters_to_fit = ones(n_parameters, 1, 'int32');
end

% now only weights and user_info could be not given (empty matrix)

%% type checks

% data, weights (if given), initial_parameters are all single
if ~isa(data, 'single')
    data = single(data);
end

if ~isempty(weights)
    if ~isa(weights, 'single')
        weights = single(weights);
    end
end

if ~isa(initial_parameters, 'single')
    initial_parameters = single(initial_parameters);
end

% parameters_to_fit is int32 (cast to int32 if incorrect type)
% if ~isa(parameters_to_fit, 'int32')
%     parameters_to_fit = int32(parameters_to_fit);
% end
if ~isa(parameters_to_fit, 'int32')
    parameters_to_fit = int32(parameters_to_fit);
end

% max_n_iterations must be int32 (cast if incorrect type)
if ~isa(max_n_iterations, 'int32')
    max_n_iterations = int32(max_n_iterations);
end

% tolerance must be single (cast if incorrect type)
if ~isa(tolerance, 'single')
    tolerance = single(tolerance);
end

% we don't check type of user_info, but we extract the size in bytes of it
if ~isempty(user_info)
    user_info_info = whos('user_info');
    user_info_size = user_info_info.bytes;
else
    user_info_size = 0;
end


%% run Gpufit taking the time
tic;
[parameters, states, chi_squares, n_iterations] ...
    = GpufitMex_1(data, weights, n_fits, n_points, tolerance, max_n_iterations, estimator_id, initial_parameters, parameters_to_fit, model_id, n_parameters, user_info, user_info_size);

time = toc;

% reshape the output parameters array to have dimensions
% (n_parameters,n_fits)
parameters = double(reshape(parameters,n_parameters,n_fits));
%convert results into double
states = double(states);
chi_squares = double(chi_squares);
n_iterations = double(n_iterations);

end
