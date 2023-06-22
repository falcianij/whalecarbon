function [n_grid] = whaleDistribution(fluxMatrix, sdmMaps, oceanMasks)
    % Computes surface species abundance distributions for whales in given years
    %
    %   Input data must be matrices of the following dimensions:
    %       - fluxMatrix[ocean, whale, year]
    %       - sdmMaps[lat, long, whale]
    %       - oceanMasks[lat, long, ocean]
    %
    %   Output is a matrix of abundance per cell of the following dimensions:
    %       - n_grid[lat, long, whale, year]
    

    n_grid = zeros(size(sdmMaps, 1), ...   % lat 
        size(sdmMaps, 2), ...              % long
        size(sdmMaps, 3), ...              % whale 
        size(fluxMatrix, 3));              % year

    for year = 1:size(n_grid, 4)
        for whale = 1:size(n_grid, 3)
            for ocean = 1:size(oceanMasks, 3)
                % abundance * (cell suitability of ocean / total suitability of ocean)
                n_cell = fluxMatrix(ocean, whale, year) .* (sum(sum(sdmMaps(:,:,whale) .* oceanMasks(:,:,ocean))) .\ (sdmMaps(:,:,whale) .* oceanMasks(:,:,ocean)) );
                n_grid(:,:, whale, year) = n_grid(:,:, whale, year) + n_cell;
            end
        end
    end

end