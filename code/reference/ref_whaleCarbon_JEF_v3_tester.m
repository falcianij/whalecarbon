%% Filepath finagling
wd = strcat(pwd, '\..\'); % working directory


%% Initialization
name_list = ["blue", "bowhead", "bryde", "fin", "gray", "humpback", "minke", "right", "sei"]; % index 1
oceans_list = ["NA", "NP", "SH", "AO", "NWP", "NEP"]; % index 2
years = ["1900", "2001"]; % index 3


%% Data finagling
ctl = load(strcat(wd, 'data_out\CTL.mat')); 
output = ctl.output; grid = output.grid;
abundanceMatrix = load(strcat(wd, 'data_out\abundanceMatrix.mat')).abundanceMatrix  ;
sdmMaps = load(strcat(wd, 'data_out\sdmMaps.mat')).sdmMaps;
oceanMasks = load(strcat(wd, 'data_out\oceanMasks.mat')).oceanMasks;


%% Mapmaking
distrib = whaleDistribution(abundanceMatrix, sdmMaps, oceanMasks);


%% Preparation of the transport matrix
grid = output.grid;
TR = output.TR; % yr^-1
msk = output.msk;
M3d = output.M3d; % land = 0, ocean = 1
VOL = grid.DXT3d .* grid.DYT3d .* grid.DZT3d; % volume, if you want area just take last term out
V = VOL(msk.pkeep);

m = size(TR,1);
sink = zeros(m,1);
sink(1:length(msk.hkeep)) = 1e10; % instantaneous surface SINK
SSINK = spdiags(sink,0,m,m);
A = TR-SSINK; % transport + sink in surface

fg = distrib(:,:,8,2)'; % transposed! % needs to be rectified s.t. it's deadfall (flux)


%% Calculations
y = -90:2:90; ny = length(y);
x = 0:2:360; nx = length(x);
z = grid.zt; nz = length(z);
[latg,long] = meshgrid(y, x);
Q = zeros();

figure(1);
clf;
ax = axesm('eqdazim','Frame','on','MapLatLimit',[-90 -30],'Origin', [-90 -180 0],'FLineWidth',0.5);
axis off; gridm off; framem on;
colormap(cool)
surfm(latg,long,fg); %Surface plot, matlab function, lat and long (y,x), fg is the flux that's being produced 
geoshow('landareas.shp','FaceColor','#EEE8AA');
cbar = colorbar; cbar.Label.String = 'Export flux [gC m^{-2} year^{-1}]';

%% Martin curve for flux attenuarion
z0 = 20; b = 0.5; flux = @(z) (z/z0).^(-b);

zw = grid.zw';
zwstar = [zw; zw(end) + grid.dzt(end)];
dz = zwstar(2:end) - zwstar(1:end-1);

fzstar = flux(zwstar);
source = -(fzstar(2:end) - fzstar(1:end-1))./dz; 
source(1) = 0;

bottomindex = sum(M3d,3); % bottom at zwstar(bottomindex + 1), bottomindex = 0 => land 
fmsk = fg>0;

figure(5); clf;
subplot(4,1,1)
imagesc(bottomindex); axis xy;
subplot(4,1,2)
imagesc(fg'); axis xy;
subplot(4,1,3)
imagesc(fmsk'); axis xy;
source = reshape(source,[1 1 nz]);

Qz = zeros(size(M3d));
Qb = zeros(size(M3d));

[latindex,lonindex] = find(fg' > 0);
for i = 1:length(lonindex)
    lati = latindex(i);
    loni = lonindex(i); 
    if loni < 180
        s = source.*M3d(lati,loni,:);
        boti = bottomindex(lati,loni);
        if boti > 0
            Qz(lati,loni,:) = fg(loni,lati)*s;
            Qz(lati,loni,boti) = 0;
            Qb(lati,loni,boti) = fg(loni,lati)*fzstar(boti)/dz(boti);
        end
    end
end
qb = sum(Qb,3);
subplot(4,1,4)
imagesc(qb); axis xy;

%%

% sall = repmat(source,[nx ny 1]); % the source form function based on Martin curve
% fall = repmat(fg,[1 1 nz]);  % the surface flux 
% Q = sall.*fall;

[lonq,latq,zq] = meshgrid(x,y,z);
Q = Qz + Qb;
i = 11; q = squeeze(Q(latindex(i),lonindex(i),:)); figure(7); clf; stairs(q,-zwstar(2:end))
q_ocim = Q(msk.pkeep); 
q_ocim(isnan(q_ocim)) = 0;
export = V'*q_ocim / 1e15 % [PgC / yr]
tic
cseq = -A \q_ocim;
toc
TotCseq = V'*cseq / 1e15 % grams to Pg
seqtime = TotCseq / export

% Re interpolate to Q
qocim = 0*M3d+NaN; % make a 3-d array of NaN
qocim(msk.pkeep) = q_ocim;
qocim(isnan(qocim)) = 0;

C_eq = 0*M3d+NaN;
C_eq(msk.pkeep) = cseq;
C_eq(isnan(C_eq)) = 0;

%% A simple plot
cc = sum(C_eq.*grid.DZT3d,3); cc = [cc, cc(:,end)]; % gC/m2



figure(2);
clf;
cc(cc<0) = 0;
subplot(2,1,1);
axesm('MapProjection','robinson','Origin',[0 270 0])
axis off; gridm off; framem on;
title('Sequestered DIC (gC/m^{-2})')
surfm(latg,long,cc');
colormap(cool);
geoshow('landareas.shp','FaceColor','#EEE8AA')
colorbar
subplot(2,1,2)
axesm('MapProjection','robinson','Origin',[0 270 0])
axis off; gridm off; framem on;
title('Surface flux (gC m^{-2} year^{-1})')
%colormap(cool)
surfm(latg,long,fg);
geoshow('landareas.shp','FaceColor','#EEE8AA')
colorbar

