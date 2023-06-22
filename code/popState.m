function [popState] = popState(abundanceMatrix, whaleParams)
    % Computes deadfall input flux for each whale species given the "steady state" abundance
    %
    %   Input data must be matrices of the following dimensions:
    %       - abundance_matrix[ocean, whale, year]
    %       - whaleParams[param, whale], where params = [max_age, mature_age, s_juvenile, s_adult, minf_male, minf_female, k_male, k_female, a0_male, a0_female]
    %
    %   Output is a matrix of flux per cell of the following dimensions:
    %       - fluxMatrix[ocean, whale, year]
    

    whaleParams(1:2, :) = round(whaleParams(1:2, :));

    nMatrix = zeros(size(abundanceMatrix));
    biomassMatrix = zeros(size(abundanceMatrix));
    fluxMatrix = zeros(size(abundanceMatrix));

    for year = 1:size(abundanceMatrix, 3)
        for whale = 1:size(abundanceMatrix, 2)
            for ocean = 1:size(abundanceMatrix, 1)
                
                % Parameters
                K = abundanceMatrix(ocean,whale,year);
                max_age = whaleParams(1, whale); mature_age = whaleParams(2, whale);
                s_juvenile = whaleParams(3, whale); s_adult = whaleParams(4, whale);
                minf_male = whaleParams(5, whale); minf_female = whaleParams(6, whale);
                k_male = whaleParams(7, whale); k_female = whaleParams(8, whale);
                % a0_male = whaleParams(9, whale); a0_female = whaleParams(10, whale); % this doesn't make sense, so I am revising. If making manuscript, rectify here!
                a0_male = 0; a0_female = 0;

                % von Bertalanffy growth functions
                m_male = minf_male .* (1 - exp(-k_male .* ((1:max_age) - a0_male) )); % male von Bertalanffy growth array
                m_female = minf_female .* (1 - exp(-k_female .* ((1:max_age) - a0_female) )); % female von Bertalanffy growth array

                % Age, biomass, and death distributions
                survivalStruct = [1, ...
                    s_juvenile .* (s_adult.^((2:max_age)-2))]; % percent of N1 that survives to j age
                n1 = K ./ (1 + sum( survivalStruct(2:max_age) )); % N1 = number of year 0 juveniles
                nStruct = n1 * survivalStruct; % number n of whales at j age
                bStruct = (0.5 .* nStruct .* m_male) + (0.5 .* nStruct .* m_female); % total biomass b of whales at j age
                
                deathStruct = [1-s_juvenile, ... % percent of N1 that dies as a juvenile
                    s_juvenile.*(1-s_adult).*(s_adult.^((2:(max_age-1))-2)), ... % percent of N1 that dies at j age
                    s_juvenile.*(s_adult.^(max_age-2))]; % percent of N1 that dies at max_age
                nDeathStruct = deathStruct .* n1; % number n of whales that die at j age
                bDeathStruct = (0.5 .* nDeathStruct .* m_male) + (0.5 .* nDeathStruct .* m_female); % total biomass b of whales that die at j age

                % Output
                nMatrix(ocean, whale, year) = sum(nStruct);
                biomassMatrix(ocean, whale, year) = sum(bStruct);
                fluxMatrix(ocean, whale, year) = sum(bDeathStruct);

            end
        end
    end
    
    popState.abundance = nMatrix;
    popState.biomass = biomassMatrix;
    popState.flux = fluxMatrix;

end