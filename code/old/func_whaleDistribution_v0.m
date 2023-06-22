%% Filepath finagling
wd = strcat(pwd, '\..\'); % working directory, dependent on file structure


%% Initialization
name_list = ["blue", "bowhead", "bryde", "fin", "gray", "humpback", "minke", "right", "sei"]; % index 1
oceans_list = ["NA", "NP", "SH", "AO", "NWP", "NEP"]; % index 2
years = ["1900", "2001"]; % index 3


%% Data finagling
load(strcat(wd, 'data_out\CTL.mat')); grid = output.grid;
abundance_matrix = load(strcat(wd, 'data_out\abundance.mat')).abundance_matrix  ;
sdm_maps = load(strcat(wd, 'data_out\sdm_maps.mat')).sdm_matrix_adj;
ocean_mask = load(strcat(wd, 'data_out\ocean_mask.mat')).ocean_mask;


%% dd
NAocean_gridArea = grid.DXT3d(:,:,1) .* grid.DYT3d(:,:,1) .* ocean_mask(:,:,1);
NAocean_totArea = sum(sum(NAocean_gridArea)); % total area North Atlantic

sdm_grid = sdm_maps .* ocean_mask(:,:,1);

abundance_per_area_BLUE_NA_1900 = abundance_matrix(2, 4, 1); % 1 = blue; 1 = NA; 1 = 1900
suitability_per_cell_BLUE_NA_1900 = sum(sum(sdm_maps(:,:,2) .* ocean_mask(:,:,4))) .\ (sdm_maps(:,:,2) .* ocean_mask(:,:,4));
abundance_per_cell_BLUE_NA_1900 = abundance_per_area_BLUE_NA_1900 .* suitability_per_cell_BLUE_NA_1900;
sum(sum(abundance_per_cell_BLUE_NA_1900))
abundance_per_area_BLUE_NA_1900

spc = sum(sum(sdm_maps(:,:,:) .* ocean_mask(:,:,3))) .\ (sdm_maps(:,:,:) .* ocean_mask(:,:,3));

%%
whaleDistribution(abundance_matrix, sdm_maps, ocean_mask);

%% Functions
function [n_grid] = whaleDistribution(abundance_matrix, sdm_maps, ocean_mask)
    % Computes species abundance distributions for whales in given years
    %
    %   Input data must be matrices of the following dimensions:
    %       - abundance_matrix[ocean, whale, year]
    %       - sdm_maps[lat, long, whale]
    %       - ocean_mask[lat, long, ocean]
    %
    %   Output is a matrix of abundance per cell of the following
    %   dimensions:
    %       - n_grid[lat, long, whale, year]
    

    n_grid = zeros(size(sdm_maps, 1), ...   % lat 
        size(sdm_maps, 2), ...              % long
        size(sdm_maps, 3), ...              % whale 
        size(abundance_matrix, 3));         % year

    for year = 1:size(n_grid, 4)
        for whale = 1:size(n_grid, 3)
            for ocean = 1:size(ocean_mask, 3)
   
                % 
                n_cell = abundance_matrix(ocean, whale, year) .* (sum(sum(sdm_maps(:,:,whale) .* ocean_mask(:,:,ocean))) .\ (sdm_maps(:,:,whale) .* ocean_mask(:,:,ocean)) );
                n_grid(:,:, whale, year) = n_grid(:,:, whale, year) + n_cell;

            end
        end
    end
end

