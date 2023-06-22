function [fluxMatrix] = fluxMatrix(abundanceMatrix, whaleParams)
    % Computes deadfall input flux for each whale species given the "steady state" abundance
    %
    %   Input data must be matrices of the following dimensions:
    %       - abundance_matrix[ocean, whale, year]
    %       - whaleParams[param, whale], where params = [max_age, mature_age, s_juvenile, s_adult, minf_male, minf_female, k_male, k_female, a0_male, a0_female]
    %
    %   Output is a matrix of flux per cell of the following dimensions:
    %       - fluxMatrix[ocean, whale, year]
    

    whaleParams(1:2, :) = round(whaleParams(1:2, :));

    fluxMatrix = zeros(size(abundanceMatrix));

    for year = 1:size(abundanceMatrix, 3)
        for whale = 1:size(abundanceMatrix, 2)
            for ocean = 1:size(abundanceMatrix, 1)
                
                max_age = whaleParams(1, whale); mature_age = whaleParams(2, whale);
                s_juvenile = whaleParams(3, whale); s_adult = whaleParams(4, whale);
                minf_male = whaleParams(5, whale); minf_female = whaleParams(6, whale);
                k_male = whaleParams(7, whale); k_female = whaleParams(8, whale);
                % a0_male = whaleParams(9, whale); a0_female = whaleParams(10, whale); % this doesn't make sense, so I am revising. If making manuscript, rectify here!
                a0_male = 0; a0_female = 0;

                syms k;
                nk = symsum(s_juvenile * (s_adult^(k-2)), k, 2, max_age); % number of spawning age
                n1 = abundanceMatrix(ocean, whale, year) ./ (1 + nk); % number of calves
                
                m_male = minf_male .* (1 - exp(-k_male .* ((1:max_age) - a0_male) )); % male von Bertalanffy growth array
                m_female = minf_female .* (1 - exp(-k_female .* ((1:max_age) - a0_female) )); % female von Bertalanffy growth array

                % calf deaths (male+female) TIMES calf mass; terminal deaths (male+female) TIMES asymptotic mass; noncalf deaths (male+female) TIMES noncalf masses
                % all in units of mass given in whaleParams (minf_male and minf_female)
                flux_calf = ( 0.5 * n1 * (1-s_juvenile) * m_male(1) ) + ( 0.5 * n1 * (1-s_juvenile) * m_female(1) );
                flux_terminal = ( 0.5 * n1 * (s_juvenile * (s_adult^(max_age-2))) * m_male(max_age) ) + ( 0.5 * n1 * (s_juvenile * (s_adult^(max_age-2))) * m_female(max_age) );
                flux_noncalf = symsum(( 0.5 * n1 * ((1-s_adult) * s_juvenile * (s_adult^(k-2))) * (minf_male .* (1 - exp(-k_male .* (k - a0_male) ))) ) + ...
                    ( 0.5 * n1 * ((1-s_adult) * s_juvenile * (s_adult^(k-2))) * (minf_female .* (1 - exp(-k_female .* (k - a0_female) ))) ), ...
                    k, 2, max_age-1);

                %flux_noncalf = 0;
                %for age = 2:(max_age-1) % extremely inefficient workaround
                    %flux_noncalf = flux_noncalf + ( 0.5 * n1 * ((1-s_adult) * s_juvenile * (s_adult^(age-2))) * m_male(age) ) + ( 0.5 * n1 * ((1-s_adult) * s_juvenile * (s_adult^(age-2))) * m_female(age) );
                %end

                fluxMatrix(ocean, whale, year) = flux_calf + flux_terminal + flux_noncalf;

            end
        end
    end

end