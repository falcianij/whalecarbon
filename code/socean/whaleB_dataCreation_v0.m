%% Filepath finagling
wd = strcat(pwd, '\..\..\'); % working directory, dependent on file structure


%% Transport matrix
load(strcat(wd, 'data\CTL.mat')); % lower resolution - less than a minute
grid = output.grid;


%% Building whale dataset
% rows are LATITUDE (0<y<180, increment 2° for 91 elements); columns are LONGITUDE (0<x<358, increment 2° for 180 elements)
globe_area_grid = grid.DXT3d(:,:,1) .* grid.DYT3d(:,:,1); % global surface area per cell
so_grid = zeros(size(globe_area_grid)); so_grid(1:24,:) = 1; so_grid = so_grid .* output.M3d(:,:,1); % SOcean = 1, everything else = 0
%so_area_grid = globe_area_grid .* so_grid; % surface area of SOcean per cell
%so_area_total = sum(sum(so_area_grid)); % total surface area of SOcean
%so_percent_grid = so_area_grid / so_area_total; % SOcean cell area weights
whaleB_so_grid = 0.201 * so_grid; % gC * m^-2 * yr^-1   %2.6 * so_percent_grid; 


%% Saving whale dataset
writematrix(whaleB_so_grid, strcat(wd, 'data\whaleB_gCm-2.csv'))