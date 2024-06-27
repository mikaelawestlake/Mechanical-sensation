function [q, lambda1_q, lambda2_q] = Newq_mult(j, k, x, y, q_og, lambda1, lambda2, shape)

%Calculates jump proportion matrix, q, for nodes surrounding removed nodes
% that are 3 cells stacked vertically (k-1,j), (k,j), (k+1,j), or
% horizontally (k, j-1), (k,j), (k, j+1)


% For verticlly stacked removed nodes
if strcmp(shape,'vert') == 1

    if k == x-2 && j == y-1 %top left corner
        q = q_og;
        q(3,3) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x-2 && j == y %top centre
        q = q_og;
        q(3,2) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x-2 && j == y+1 %top right corner
        q = q_og;
        q(3,1) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x-1 && j == y+1 %upper right side
        q = q_og;
        q(2,1) = 0;
        q(3,1) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x && j == y+1 % centre right 
        q = q_og;
        q(2,1) = 0;
        q(1,1) = 0;
        q(3,1) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x+1 && j == y+1 %lower right 
        q = q_og;
        q(1,1) = 0;
        q(2,1) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x+2 && j == y+1 %right bottom corner
        q = q_og;
        q(1,1) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x+2 && j == y %bottom centre
        q = q_og;
        q(1,2) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x+2 && j == y-1  %left bottom corner 
        q = q_og;
        q(1,3) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x+1 && j == y-1 %left lower side
        q = q_og;
        q(2,3) = 0;
        q(1,3) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x && j == y-1 %left centre side
        q = q_og;
        q(1,3) = 0;
        q(2,3) = 0;
        q(3,3) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;  
    
    elseif k == x-1 && j == y-1 %left upper side
        q = q_og;
        q(2,3) = 0;
        q(3,3) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x-1 && j == y %removed node (top)
        q = q_og;
        q(2,2) = 1; 
        lambda1_q = 0;
        lambda2_q = 0;
    
    elseif k == x && j == y %removed node (centre)
        q = q_og;
        q(2,2) = 1;
        lambda1_q = 0;
        lambda2_q = 0;
    
    elseif k == x+1 && j == y %removed node (bottom)
        q = q_og; 
        q(2,2) = 1; 
        lambda1_q = 0;
        lambda2_q = 0;

    else 
    q=q_og; 
    lambda1_q = lambda1;
    lambda2_q = lambda2;
  
    end

% For horizontally stacked removed nodes
elseif strcmp(shape,'horz')

    if k == x-1 && j == y-2 %top left corner
        q = q_og;
        q(3,3) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x-1 && j == y-1 %top left centre 
        q = q_og;
        q(3,2) = 0;
        q(3,3) = 0; 
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x-1 && j == y %top centre
        q = q_og;
        q(3,1) = 0;
        q(3,2) = 0;
        q(3,3) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x-1 && j == y+1 % top right centre 
        q = q_og;
        q(3,2) = 0;
        q(3,1) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x-1 && j == y+2 % top right corner
        q = q_og;
        q(3,1) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x && j == y+2 %centre right 
        q = q_og;
        q(2,1) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x+1 && j == y+2 %right bottom corner
        q = q_og;
        q(1,1) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x+1 && j == y+1 %bottom right centre 
        q = q_og;
        q(1,2) = 0;
        q(1,1) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x+1 && j == y  %bottom centre
        q = q_og;
        q(1,3) = 0;
        q(1,2) = 0;
        q(1,1) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x+1 && j == y-1 %left lower centre 
        q = q_og;
        q(1,2) = 0;
        q(1,3) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;
    
    elseif k == x+1 && j == y-2 %left lower 
        q = q_og;
        q(1,3) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;  
    
    elseif k == x && j == y-2 %left side 
        q = q_og;
        q(2,3) = 0;
        lambda1_q = lambda1;
        lambda2_q = lambda2;

    elseif k == x && j == y-1 % removed node (left)
        q = q_og;
        q(2,2) = 1;
        lambda1_q = 0;
        lambda2_q = 0;
    
    elseif k == x && j == y %removed node (centre)
        q = q_og;
        q(2,2) = 1;
        lambda1_q = 0;
        lambda2_q = 0;
    
    elseif k == x && j == y+1 %removed node (right)
        q = q_og; 
        q(2,2) = 1; 
        lambda1_q = 0;
        lambda2_q = 0;

        else 
        q=q_og; 
        lambda1_q = lambda1;
        lambda2_q = lambda2;
  
    end 

end

if sum(q, 'all') ~= 1
    q(2,2) = q(2,2) + (1 - sum(q, "all")); 
end 

end 