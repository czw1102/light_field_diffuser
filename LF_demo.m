clear;
close all;

%% Pre-processing parameters  

data_path = 'Image/';
pattern_name = 'pattern-6.tif';
image_name = 'image-1.tif';

% The desired spatial sampling
spatial_sampling_s = 256;
spatial_sampling_t = 256;

% The sub-image interval to obtain the desired angular sampling
% An approximate value dependent on the spatio-angular sampling is required
subimage_interval = 10*3;
edge_skip = 5;

% A value (usually in the range of [-0.1, 0.1]) is used to adjust the size of the support.
scaling_factor = 0;

%% Optimization parameters

% Solver parameters
solverSettings.tau = 0.001;    %sparsity parameter for TV
solverSettings.tau_n = 0.001;     %sparsity param for native sparsity
solverSettings.mu1 = 1;    %Initialize ADMM tuning params. If autotune is on, these will change
solverSettings.mu2 = 1;
solverSettings.mu3 = 1;
  
% if set to 1, auto-find mu1, mu2, mu3 every step. If set to 0, use user defined values. If set to N>1, tune for N steps then stop.
solverSettings.autotune = 1;    % default: 1
solverSettings.mu_inc = 1.1; %1.1;  % 
solverSettings.mu_dec = 1.1; %1.1;  %Inrement and decrement values for mu during autotune. Turn to 1 to have no tuning.
solverSettings.resid_tol = 1.5;   % Primal/dual gap tolerance. Lower means more frequent tuning
solverSettings.maxIter = 500; % Maximum iteration count  Default: 200
solverSettings.regularizer = 'tv';   %'TV' for 3D TV, 'native' for native. Default: TV
 
%Figures and user info
solverSettings.disp_percentile = 100;   %Percentile of max to set image scaling
solverSettings.disp_func = @(x)x;  %Function handle to modify image before display. No change to data, just for display purposes
solverSettings.disp_figs = 1;   %If set to 0, never display. If set to N>=1, show every N.
solverSettings.print_interval = 1;  %Print cost every N iterations. Default 1. If set to 0, don't print.
solverSettings.fig_num = 1;   %Figure number to display in

%% Sub-image set

pattern_base = double(imread([data_path,pattern_name]));
pattern_base = imresize(pattern_base,[spatial_sampling_t, spatial_sampling_s],'box');

image_obj = double(imread([data_path,image_name]));
image_obj = imresize(image_obj,[spatial_sampling_t, spatial_sampling_s],'box');
image_obj = mat2gray(image_obj);

% The support (encoding kernel) of the pattern
% Left and right edges of the support
pattern_base_sum = sum(pattern_base,1);
pattern_base_sum_mean = mean(pattern_base_sum);
pattern_base_sum_ind = find(pattern_base_sum>pattern_base_sum_mean/10);
edge_left_round = pattern_base_sum_ind(1)-edge_skip;
edge_right_round = pattern_base_sum_ind(end)+edge_skip;

edge_scaling_size = round(spatial_sampling_s*scaling_factor);
edge_left_round_scaling = edge_left_round-edge_scaling_size;
edge_right_round_scaling = edge_right_round+edge_scaling_size;

edge_lr_interval = edge_right_round_scaling-edge_left_round_scaling+1;
angular_sampling_u = edge_lr_interval/subimage_interval;

if mod(floor(angular_sampling_u),2)
    edge_lr_interval_diff = (floor(angular_sampling_u)+1)*subimage_interval-edge_lr_interval;
    edge_left = edge_left_round_scaling-ceil((edge_lr_interval_diff-1)/2);
    edge_right = edge_right_round_scaling+floor((edge_lr_interval_diff+1)/2);
else
    edge_lr_interval_diff = edge_lr_interval-floor(angular_sampling_u)*subimage_interval;
    edge_left = edge_left_round_scaling+ceil((edge_lr_interval_diff-1)/2);
    edge_right = edge_right_round_scaling-floor((edge_lr_interval_diff+1)/2);
end

% Up and down edges of the support
pattern_base_sum = sum(pattern_base,2);
pattern_base_sum_mean = mean(pattern_base_sum);
pattern_base_sum_ind = find(pattern_base_sum>pattern_base_sum_mean/10);
edge_up_round = pattern_base_sum_ind(1)-edge_skip;
edge_down_round = pattern_base_sum_ind(end)+edge_skip;

edge_scaling_size = round(spatial_sampling_t*scaling_factor);
edge_up_round_scaling = edge_up_round-edge_scaling_size;
edge_down_round_scaling = edge_down_round+edge_scaling_size;

edge_ud_interval = edge_down_round_scaling-edge_up_round_scaling+1;
angular_sampling_v = edge_ud_interval/subimage_interval;

if mod(floor(angular_sampling_v),2)
    edge_ud_interval_diff = (floor(angular_sampling_v)+1)*subimage_interval-edge_ud_interval;
    edge_up = edge_up_round_scaling-ceil((edge_ud_interval_diff-1)/2);
    edge_down = edge_down_round_scaling+floor((edge_ud_interval_diff+1)/2);
else
    edge_ud_interval_diff = edge_ud_interval-floor(angular_sampling_v)*subimage_interval;
    edge_up = edge_up_round_scaling+ceil((edge_ud_interval_diff-1)/2);
    edge_down = edge_down_round_scaling-floor((edge_ud_interval_diff+1)/2);
end

% The desired angular sampling
angular_sampling_u = (edge_right-edge_left+1)/subimage_interval;
angular_sampling_v = (edge_down-edge_up+1)/subimage_interval;

% The scaled pattern
pattern_base_new = zeros(size(pattern_base));
pattern_base_new(edge_up_round_scaling:edge_down_round_scaling, edge_left_round_scaling:edge_right_round_scaling) = ...
    imresize(pattern_base(edge_up_round:edge_down_round, edge_left_round:edge_right_round), [edge_ud_interval, edge_lr_interval], 'box');

% The sub-image set according to the spatio-angular sampling
subimage_set = zeros(spatial_sampling_t,spatial_sampling_s,angular_sampling_v,angular_sampling_u);
ind = 0;
for j = 1:angular_sampling_u
    for i = 1:angular_sampling_v
        ind = ind+1;
        subimage_set(edge_up-1+(i-1)*subimage_interval+(1:subimage_interval),edge_left-1+(j-1)*subimage_interval+(1:subimage_interval),i,j)...
            = pattern_base(edge_up-1+(i-1)*subimage_interval+(1:subimage_interval),edge_left-1+(j-1)*subimage_interval+(1:subimage_interval));
    end
end

%% Reconstruction

[lf_obj,opt_param] = LF_ADMM_solver(subimage_set,image_obj,solverSettings);

save_file = [pattern_name(1:end-4),'-',image_name(1:end-4),'-',...
    num2str(angular_sampling_v),'x',num2str(angular_sampling_u),'x',num2str(spatial_sampling_t),'x',num2str(spatial_sampling_s),'.mat'];
save([data_path,save_file],'lf_obj','opt_param');
