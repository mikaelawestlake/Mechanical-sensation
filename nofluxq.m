function [q] = nofluxq(k, j, q_og, K, M)

% Calculates the jump proporiton matrix, q, based on current location in
% lattice for node along the no flux boundary's (Upper, left, lower)

% UPPER LEFT CORNER 
if k == 2 && j == 2

    q = [ 0 0 0; 0 0.75 0.1; 0 0.1 0.05]; 

% % LEFT BOUNDARY - SIDES
elseif 2 <= k && k < K-1  && j == 2

    q = [0 0.1 0.05; 0 0.65 0.05; 0 0.1 0.05];

% LOWER LEFT CORNER 
elseif k == K-1 && j == 2

    q = [0 0.1 0.05; 0 0.75 0.1; 0 0 0]; 
    
% LOWER BOUNDARY NO FLUX - SIDES
elseif k == K-1 && j ~= 2

    q = [0.05 0.05 0.05; 0.1 0.65 0.1 ; 0 0 0];


% UPPER BOUNDARY NO FLUX - SIDES
elseif k == 2 && j ~=2

    q = [0 0 0; 0.1 0.65 0.1; 0.05 0.05 0.05];

else
    q = q_og;
end 

end 