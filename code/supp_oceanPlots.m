%% Filepath finagling
wd = strcat(pwd, '\..\'); % working directory
%codePath = addpath(strcat(wd, 'code'), '-end');

%% Initialization
name_list = ["blue", "bowhead", "bryde", "fin", "gray", "humpback", "minke", "right", "sei"]; % index 1
oceans_list = ["NA", "NP", "SH", "AO", "NWP", "NEP"]; % index 2
years = ["1900", "2001"]; % index 3


%% Data finagling
CTL = load(strcat(wd, 'data_out\CTL.mat')); output = CTL.output; grid = output.grid;
abundanceMatrix = load(strcat(wd, 'data_out\abundanceMatrix.mat')).abundanceMatrix;
sdmMaps = load(strcat(wd, 'data_out\sdmMaps.mat')).sdmMaps;
oceanMasks = load(strcat(wd, 'data_out\oceanMasks.mat')).oceanMasks;
whaleParams = load(strcat(wd, 'data_out\whaleParams.mat')).whaleParams;

% Predictor raster of interest
[A,R] = readgeoraster(strcat(wd, 'data\gis\clipped\chl_clipped_res.tif'));
%info = geotiffinfo(strcat(wd, 'data\gis\clipped\vel_clipped_res.tif'));
[lat,lon] = geographicGrid(R);


%% IWC data plotting
y = -90:2:90; ny = length(y);
x = 0:2:360; nx = length(x);
z = grid.zt; nz = length(z);
[latg,long] = meshgrid(y, x);

figure(1);
clf;
axesm('MapProjection','robinson','Origin',[0 270 0])
axis off; gridm off; framem on;
title('Chlorophyll concentration (mg·m^{-3})') %(mg·m^{-3})
surfm(lat, lon, log(A)); % for predictor raster of interest
%surfm(latg,long,oceans_grid');
colormap(cool);
geoshow('landareas.shp','FaceColor','#EEE8AA')
colorbar


%% Predictor variable plotting
%y = -90:2:90; ny = length(y);
%x = 0:2:360; nx = length(x);
%z = grid.zt; nz = length(z);
%[latg,long] = meshgrid(y, x);

%figure(1);
%clf;
%axesm('MapProjection','robinson','Origin',[0 270 0])
%axis off; gridm off; framem on;
%title('Chlorophyll concentration (mg·m^{-3})') %(mg·m^{-3})
%surfm(lat, lon, log(A)); % for predictor raster of interest
%surfm(latg,long,oceans_grid');
%colormap(cool);
%geoshow('landareas.shp','FaceColor','#EEE8AA')
%colorbar


%% Ocean delineation
oceanMasks(:,:,2) = oceanMasks(:,:,2) * 0;
oceanMasks(:, :, 3) = oceanMasks(:,:,3) * 2;
oceanMasks(:, :, 4) = oceanMasks(:,:,4) * 3;
oceanMasks(:, :, 5) = oceanMasks(:,:,5) * 4;
oceanMasks(:, :, 6) = oceanMasks(:,:,6) * 5;

oceans_grid = sum(oceanMasks, 3);


%%
%y = -90:2:90; ny = length(y);
%x = 0:2:360; nx = length(x);
%z = grid.zt; nz = length(z);
%[latg,long] = meshgrid(y, x);

%figure(1);
%clf;
%axesm('MapProjection','robinson','Origin',[0 270 0])
%axis off; gridm off; framem on;
%title('Sequestered DIC (gC/m^{-2})')
%surfm(latg,long,oceans_grid');
%colormap(cool);
%geoshow('landareas.shp','FaceColor','#EEE8AA')
%colorbar