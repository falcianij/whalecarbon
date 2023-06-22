function [cState] = transport(cparams, CTL, fdistribMap)
    % Computes the 3D map of carbon sequestration given an input flux map and sequestration paramters
    %
    %   Input data must of the following dimensions:
    %       - params[parameters] = [z0, b]
    %       - CTL matrix (provided by DeVries et al.)
    %       - fdistribMap[lat, long]
    %
    %   Output is a struct of of the following:
    %       - State.export = carbon exported (PgC / yr)
    %       - cState.totCseq = total carbon sequestered (Pg)
    %       - cState.seqtime = residence time of sequestered carbon (yr)
    %       - cState.cGrid = carbon remineralization profile [lat, long, depth];
    %       - cState.cBottom = carbon accumulated at the seafloor [lat, long];
    

    % Preparation of the transport matrix
    output = CTL.output;
    grid = output.grid;
    TR = output.TR; % yr^-1     % transport
    msk = output.msk;
    M3d = output.M3d; % land = 0, ocean = 1
    VOL = grid.DXT3d .* grid.DYT3d .* grid.DZT3d; % volume of every cell
    V = VOL(msk.pkeep); % volume of oceans only

    m = size(TR,1);
    sink = zeros(m,1);
    sink(1:length(msk.hkeep)) = 1e10; % instantaneous surface SINK
    SSINK = spdiags(sink,0,m,m);
    A = TR-SSINK; % transport + sink in surface

    fg = fdistribMap'; % transposed!
    

    % Calculations
    y = -90:2:90; ny = length(y); % dimensions of lat (for everything; from CTL)
    x = 0:2:360; nx = length(x); % dimensions of long (for everything; from CTL)
    z = grid.zt; nz = length(z); % dimensions of depth
    Q = zeros();


    % Martin curve for flux attenuarion
    z0 = cparams(1); b = cparams(2);
    flux = @(z) (z/z0).^(-b);

    zw = grid.zw';
    zwstar = [zw; zw(end) + grid.dzt(end)];
    dz = zwstar(2:end) - zwstar(1:end-1);

    fzstar = flux(zwstar);
    source = -(fzstar(2:end) - fzstar(1:end-1))./dz; 
    source(1) = 0;
    source = reshape(source,[1 1 nz]);
    
    bottomindex = sum(M3d,3); % bottom at zwstar(bottomindex + 1), bottomindex = 0 => land 


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


    Q = Qz + Qb;
    q_ocim = Q(msk.pkeep);
    q_ocim(isnan(q_ocim)) = 0;
    export = V'*q_ocim / 1e15; % [PgC / yr]
    %tic
    cseq = -A \ q_ocim;
    %toc
    totCseq = V'*cseq / 1e15; % grams to Pg
    seqtime = totCseq / export;

    C_eq = 0*M3d+NaN;
    C_eq(msk.pkeep) = cseq;
    C_eq(isnan(C_eq)) = 0;

    % Bottom C
    cc = sum(C_eq.*grid.DZT3d,3); cc = [cc, cc(:,end)]; % gC/m2

    
    % Outputs
    cState.export = export; % [PgC / yr]
    cState.totCseq = totCseq; % grams to Pg
    cState.seqtime = seqtime; % yr
    cState.cGrid = C_eq;
    cState.cBottom = cc;
end