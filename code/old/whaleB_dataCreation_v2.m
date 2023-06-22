%% Filepath finagling
wd = strcat(pwd, '\..\'); % working directory, dependent on file structure


%% Transport matrix
load(strcat(wd, 'data\CTL.mat')); % lower resolution - less than a minute
image = Tiff(strcat(wd, 'data\obis\right\Right_sdm_matlab.tif'));
whaleB_grid = read(image);
whaleB_grid(whaleB_grid <= -9999) = 0;

whaleB_grid_adj = zeros(size(whaleB_grid));
whaleB_grid_adj(:, 1:90) = whaleB_grid(:, 91:180);
whaleB_grid_adj(:, 91:180) = whaleB_grid(:, 1:90);

whaleB_grid_adj = flip(whaleB_grid_adj);


%% Building whale dataset
% rows are LATITUDE (0<y<180, increment 2° for 91 elements); columns are LONGITUDE (0<x<358, increment 2° for 180 elements)
whaleB_grid_adj = whaleB_grid_adj .* output.M3d(:,:,1); % masking input data to land, water defined by transport matrix


%% Saving whale dataset
writematrix(whaleB_grid_adj, strcat(wd, 'data\right_whaleB_gCm-2.csv'))