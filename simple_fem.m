
function [SED] = simple_fem(nelx, nely, forces, applied, psi_bar)

% SIMPLE 2D Linear Elastic Finite Element method code
% Vivien Challis, 2024
% Adapted from supplementary code published 2010 and available at:
%     dx.doi.org/10.1007/s00158-009-0430-0
% The above code was adapted from the earlier code available at:
%     dx.doi.org/10.1007/s001580050176
% Follow this link to find a description of the node numbering. Both nodes
% and elements are numbered column wise from left to right.


struc = ones(nely,nelx); % ADJUST FOR THE PROBLEM
% struc(3:4,4:7)=0;
SED = zeros(nely,nelx);
[KE,~,~] = materialInfo();

 % FE-analysis
 [U, ~] = FE(struc,KE, applied);
 for ely = 1:nely
  for elx = 1:nelx
   n1 = (nely+1)*(elx-1)+ely;
   n2 = (nely+1)* elx   +ely;
   Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
   SED(ely,elx) = 0.5*max(struc(ely,elx),0.0001)*Ue'*KE*Ue;
  end
 end


 % Visualisation 
 cmap = colormap(parula(20));
surf(SED/psi_bar)
view(2)
cb = colorbar
clim([0 10])
set(gca,'colorscale','linear
axis equal
xlim([1 11])
ylim([1 31])
title('Strain Energy Density', 'FontSize', 12, 'FontName', 'times')
xlabel('Position in osteon in j (nodes)' , 'FontSize', 12, 'FontName', 'times')
ylabel('Position in osteon in k (nodes)', 'FontSize', 12, 'FontName', 'times')
set(gca, 'YDir','reverse')
yticks([1 6 11 16 21 26 31] )
% print(gcf, sprintf('SED horz out- %f .png', applied),'-dpng','-r300')


%%---- FINITE ELEMENT ANALYSIS ----
function [U, F] = FE(struc,KE, ~)
[nely,nelx] = size(struc);
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
F = sparse(2*(nely+1)*(nelx+1),1); 
U = zeros(2*(nely+1)*(nelx+1),1);
for elx = 1:nelx
 for ely = 1:nely
  n1 = (nely+1)*(elx-1)+ely;
  n2 = (nely+1)* elx   +ely;
  edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
  K(edof,edof) = K(edof,edof) + max(struc(ely,elx),0.0001)*KE;
 end
end

% Define loads and supports - ADJUST FOR THE PROBLEM:
% Boundary conditions represent the solid sitting on a hard surface.
% This line implements a spatially homogeneous constant force on the top 
% boundary:


fixeddofs = [2*nely+1];

for nx = 1:nelx
    F(2*(nx-1)*(nely+1)+2,1) = F(2*(nx-1)*(nely+1)+2,1) - 0.5*forces(nx);
    F(2*(nx)*(nely+1)+2,1) = F(2*(nx)*(nely+1)+2,1) - 0.5*forces(nx);
    fixeddofs = [fixeddofs, 2*nx*(nely+1)];
end
fixeddofs = [fixeddofs, 2*(nelx+1)*(nely+1)];

% Solving
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
U(freedofs,:) = K(freedofs,freedofs) \ F(freedofs,:);

end

%%---- MATERIAL INFORMATION ----
    function [KE,lambda,mu] = materialInfo()
% Set material parameters, find Lame values
E = 28e+09; nu = 0.3; % ADJUST FOR THE PROBLEM E in GPa, nu - poisson
lambda = E*nu/((1+nu)*(1-nu));

mu = E/(2*(1+nu)); %Lame value

% Find stiffness matrix "KE"
k =[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ...
    -1/4+nu/12 -1/8-nu/8       nu/6  1/8-3*nu/8];
KE = E/(1-nu^2)*stiffnessMatrix(k);
end

%%---- ELEMENT STIFFNESS MATRIX ----
function [K] = stiffnessMatrix(k)
% Forms stiffness matrix from first row

K=[k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
   k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
   k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
   k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
   k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
   k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
   k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
   k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
end
end