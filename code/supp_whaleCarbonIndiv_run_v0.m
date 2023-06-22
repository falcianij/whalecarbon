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

whale = 7;


year = 1;
    cState1 = transport(cparams, CTL, fdistrib(:, :, whale, year));

    cc1 = cState1.cBottom;
    fg1 = fdistrib(:, :, whale, year);

    Cseq1 = cState1.totCseq;
    totCseq1 = cState1.totCseq;
    totExport1 = cState1.export
    totSeqTime1 = totCseq1 / totExport1;

year = 2;
    cState2 = transport(cparams, CTL, fdistrib(:, :, whale, year));

    cc2 = cState2.cBottom;
    fg2 = fdistrib(:, :, whale, year);

    Cseq2 = cState2.totCseq;
    totCseq2 = cState2.totCseq;
    totExport2 = cState2.export;
    totSeqTime2 = totCseq2 / totExport2;


%cc = cState.cBottom;
%fg = fdistrib(:, :, whale, year)';
fg1 = fg1';
fg2 = fg2';
%% Plotting


y = -90:2:90; ny = length(y);
x = 0:2:360; nx = length(x);
z = grid.zt; nz = length(z);
[latg,long] = meshgrid(y, x);

minColorLimit = (min(min(fg1)));
maxColorLimit = (max(max(fg1)));

figure(2);
clf;
cc1(cc1<0) = 0;
cc2(cc2<0) = 0;
subplot(1,2,1);
axesm('MapProjection','robinson','Origin',[0 270 0])
axis off; gridm off; framem on;
%title('Sequestered DIC (gC/m^{-2})')
surfm(latg,long,fg1);
geoshow('landareas.shp','FaceColor','#EEE8AA')
%colorbar
subplot(1,2,2)
axesm('MapProjection','robinson','Origin',[0 270 0])
axis off; gridm off; framem on;
%title('Surface flux (gC m^{-2} year^{-1})')
%colormap(cool)
surfm(latg,long,fg2);
geoshow('landareas.shp','FaceColor','#EEE8AA')
%colorbar

colormap(cool);
colorbar
caxis([minColorLimit,maxColorLimit]);


figure(3);
clf;
cc1(cc1<0) = 0;
cc2(cc2<0) = 0;
subplot(1,2,1);
axesm('MapProjection','robinson','Origin',[0 270 0])
axis off; gridm off; framem on;
%title('Sequestered DIC (gC/m^{-2})')
surfm(latg,long,fg1);
geoshow('landareas.shp','FaceColor','#EEE8AA')
%colorbar
subplot(1,2,2)
axesm('MapProjection','robinson','Origin',[0 270 0])
axis off; gridm off; framem on;
%title('Surface flux (gC m^{-2} year^{-1})')
%colormap(cool)
surfm(latg,long,fg2);
geoshow('landareas.shp','FaceColor','#EEE8AA')
%colorbar

colormap(cool);
caxis([minColorLimit,maxColorLimit]);