load('CTL.mat') % lower resolution - less than a minute
%Load the CTL file which includes 3
firebrick = '#B22222';
pale_golden_rod = '#EEE8AA';
%% Preparation of the transport matrix
grid = output.grid;
TR = output.TR; % yr^-1
msk = output.msk;
M3d = output.M3d; % land = 0. ocean = 1 %3D matrix, there will be 1's at the bottom 
VOL = grid.DXT3d.*grid.DYT3d.*grid.DZT3d; %Volume, if you want area just take last term out
V = VOL(msk.pkeep);

m = size(TR,1);
sink = zeros(m,1);
sink(1:length(msk.hkeep)) = 1e10; % instantaneous Surface SINK
SSINK = spdiags(sink,0,m,m);
A = TR-SSINK; % transport + sink in surface

%% Interpolate krill flux
%Reorganizing all of the data so that it can be interpolated.. The model
%needs the neighbor data, otherwise there's a hole at 0 degrees that is
%needed... Does some cyclical math to reorient the math?
% D = readmatrix('krill_fpp_gCm-2yr-1.csv'); %reads data from the spreadsheet
% lato = D(:,1); lono = D(:,2); fo = D(:,3); %Turning data from the matrix into variables
% %Model needs to go from 0-360, right now it goes from -180 to 180
% kn = lono < 0;  %Getting the lon points that are negative
% lono(kn) = lono(kn) + 360; %Making the negative longiudes positive 
% kw = lono < 10; lonow = lono(kw)+360; %Less than 10, make a copy and put on the other side
% latow = lato(kw); fow = fo(kw);
% ke = lono > 350; lonoe = lono(ke)-360; %Have to make it overlap for the purpose of interpolation
% latoe = lato(ke); foe = fo(ke);

%Re: whale: Could also make an excel spreadsheet that has the same lat and long as the
%krill data and fill in the data regarding the cell

%Krill data has gc/m2/year
%Biomass density - total ocean size and the biomass
% B/(A T) = flux = gC/m2/year
% Area matrix is in the variable called "grid" 'grid.Areac' 
% (Area for each degree cell)
% Multiply grid.Areac by the surface layer
saw=grid.DXT3d(:,:,1).*grid.DYT3d(:,:,1).*M3d(:,:,1); %This gives you the area of each cell on the surface, only water
%so = saw(1:25,:); %Southern ocean area
sa2 = sum(saw); %Sum of the rows
satot = sum(sa2); %Total surface area (Only of water)


%how much is whale biomass in southern ocean?
sa=grid.DXT3d(:,:,1).*grid.DYT3d(:,:,1).*M3d(:,:,1); %rows are latitude and columns are longitude 
%M3d covers if you hav eland or water (0 for land)
%matrix for entire ocean
southernoceanmatrix=zeros(91,180); %building a matrix of zeros
southernoceanmatrix(1:25,:)=1; %Putting 1's at the location of the southern ocean

ocean = sa.*southernoceanmatrix; %The entire ocean, where the data is populated only at the southern ocean. 

%but this is global
sumofrows=sum(ocean);
total_area=sum(sumofrows);
percentage= ocean/total_area;
biomass=2.6*percentage; %for whole ocean

lato = [lato;latoe;latow];
lono = [lono;lonoe;lonow];
fo = [fo;foe;fow];
y = -90:2:90; ny = length(y);
x = 0:2:360; nx = length(x);
z = grid.zt; nz = length(z);
[latg,long] = meshgrid(y, x);
Q = zeros();

%fg = griddata(lato,lono,fo,latg,long);
fg = biomass';
fgNaN = isnan(fg);
fg(fgNaN)=0;
[rg,cg] = find(fg>0);



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
        [i, lati,loni,boti]
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

