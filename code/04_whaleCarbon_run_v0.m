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


%% Mapmaking
%c_percent = 0.2;
%popState = popState(abundanceMatrix, whaleParams);
%fluxMatrix = popState.flux * c_percent * 1e6; % metric tons ww -> metric tons carbon -> g carbon
%fdistrib = whaleDistribution(fluxMatrix, sdmMaps, CTL, oceanMasks); % remember: [lat, long, whale, year]

%z0 = 20; b = 0.13;
%cparams = [z0, b];

%whale = 4;
%year = 1;
%biomass = sum(popState.biomass(:, whale, year))
%flux_megatons = sum(popState.flux(:, whale, year) .* 0.5 .* c_percent ./ 1e6) % Dufort estimate
%cState = transport(cparams, CTL, fdistrib(:, :, whale, year));
%totCseq = cState.totCseq


%% Data creation
c_percent = 0.2;
z0 = 20; b = 0.13;
cparams = [z0, b];

popState = popState(abundanceMatrix, whaleParams);
fluxMatrix = popState.flux * c_percent * 1e6; % metric tons ww -> metric tons carbon -> g carbon
fdistrib = whaleDistribution(fluxMatrix, sdmMaps, CTL, oceanMasks); % remember: [lat, long, whale, year]

year = 1;
totExport = 0;
totCseq = 0;
Cseq = zeros(size(name_list));
cc = zeros(91, 181); % hard coded; dim = 2 has an extra element
fg = zeros(size(grid.DZT3d(:, :, 1)));

for whale = 1%:size(abundanceMatrix, 2)
    cState = transport(cparams, CTL, fdistrib(:, :, whale, year));

    cc = cc+ cState.cBottom;
    fg = fg + fdistrib(:, :, whale, year);

    Cseq(whale) = cState.totCseq;
    totCseq = totCseq + cState.totCseq;
    totExport = totExport + cState.export;
end

Cseq
totCseq
totExport
totSeqTime = totCseq / totExport

fg = fg';

%% Plotting
%cc = cState.cBottom;
%fg = fdistrib(:, :, whale, year)';

y = -90:2:90; ny = length(y);
x = 0:2:360; nx = length(x);
z = grid.zt; nz = length(z);
[latg,long] = meshgrid(y, x);

minColorLimit = (min(min(cc)));
maxColorLimit = (max(max(cc)));

figure(2);
clf;
cc(cc<0) = 0;
subplot(2,1,1);
axesm('MapProjection','robinson','Origin',[0 270 0])
axis off; gridm off; framem on;
%title('Sequestered DIC (gC/m^{-2})')
surfm(latg,long,cc');
colormap(cool);
geoshow('landareas.shp','FaceColor','#EEE8AA')
%caxis([0,2.17]);
colorbar
subplot(2,1,2)
axesm('MapProjection','robinson','Origin',[0 270 0])
axis off; gridm off; framem on;
%title('Surface flux (gC m^{-2} year^{-1})')
%colormap(cool)
surfm(latg,long,fg);
geoshow('landareas.shp','FaceColor','#EEE8AA')
colorbar