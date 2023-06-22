%% Filepath finagling
wd = strcat(pwd, '\..\'); % working directory


% Initialization
name_list = ["blue", "bowhead", "bryde", "fin", "gray", "humpback", "minke", "right", "sei"];
oceans_list = ["NA", "NP", "SH", "AO", "NWP", "NEP"];
years = ["1900", "2001"];


%% Assembling abundance matrix
abundanceMatrix = zeros(length(oceans_list), length(name_list), length(years));
table = readtable(strcat(wd,'data\abundanceMatrix\whale_presence_1900.csv')); % cols = name_list; rows = oceans_list
abundanceMatrix(:, :, 1) = table2array(table(:, 2:size(table, 2)));
table = readtable(strcat(wd,'data\abundanceMatrix\whale_presence_2001.csv')); % cols = name_list; rows = oceans_list
abundanceMatrix(:, :, 2) = table2array(table(:, 2:size(table, 2)));

save(strcat(wd, 'data_out\abundanceMatrix.mat'), "abundanceMatrix");


%% Assembling species distribution matrix (SDM)
sdmMatrix = zeros(91, 180, length(name_list));

sdmMatrix(:, :, 1) = read(Tiff(strcat(wd, 'data\obis\blue\blue_sdm_matlab.tif')));
sdmMatrix(:, :, 2) = read(Tiff(strcat(wd, 'data\obis\bowhead\bowhead_sdm_matlab.tif')));
sdmMatrix(:, :, 3) = read(Tiff(strcat(wd, 'data\obis\bryde\bryde_sdm_matlab.tif')));
sdmMatrix(:, :, 4) = read(Tiff(strcat(wd, 'data\obis\fin\fin_sdm_matlab.tif')));
sdmMatrix(:, :, 5) = read(Tiff(strcat(wd, 'data\obis\gray\gray_sdm_matlab.tif')));
sdmMatrix(:, :, 6) = read(Tiff(strcat(wd, 'data\obis\humpback\humpback_sdm_matlab.tif')));
sdmMatrix(:, :, 7) = read(Tiff(strcat(wd, 'data\obis\minke\minke_sdm_matlab.tif')));
sdmMatrix(:, :, 8) = read(Tiff(strcat(wd, 'data\obis\right\right_sdm_matlab.tif')));
sdmMatrix(:, :, 9) = read(Tiff(strcat(wd, 'data\obis\sei\sei_sdm_matlab.tif')));

sdmMatrix(sdmMatrix <= -9999) = 0;

sdmMaps = zeros(size(sdmMatrix));
sdmMaps(:, 1:90, :) = sdmMatrix(:, 91:180, :);
sdmMaps(:, 91:180, :) = sdmMatrix(:, 1:90, :);

sdmMaps = flip(sdmMaps);

save(strcat(wd, 'data_out\sdmMaps.mat'), "sdmMaps");


%% Assembling whale population parameters
table = readtable(strcat(wd,'data\whaleParams.csv')); % cols = name_list; rows = oceans_list
whaleParams = table2array(table(:, 2:size(table, 2)));

save(strcat(wd, 'data_out\whaleParams\whaleParams.mat'), "whaleParams");


%% Assembling ocean mask struct
oceanMasks = zeros(91, 180, length(oceans_list));

oceanMasks(:, :, 1) = flip(load(strcat(wd, 'data\oceanMasks\NA.csv')));
oceanMasks(:, :, 2) = flip(load(strcat(wd, 'data\oceanMasks\NP.csv')));
oceanMasks(:, :, 3) = flip(load(strcat(wd, 'data\oceanMasks\SH.csv')));
oceanMasks(:, :, 4) = flip(load(strcat(wd, 'data\oceanMasks\AO.csv')));
oceanMasks(:, :, 5) = flip(load(strcat(wd, 'data\oceanMasks\NWP.csv')));
oceanMasks(:, :, 6) = flip(load(strcat(wd, 'data\oceanMasks\NEP.csv')));

save(strcat(wd, 'data_out\oceanMasks.mat'), "oceanMasks");
