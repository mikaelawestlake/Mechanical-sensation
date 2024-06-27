
% Outputs at each timestep interval are saved and ca be extracted in other
% scripts for visualisaiton 

%% Establish file data directory 

if exist("datadir")==7
    rmdir("datadir", 's')
end

clear all
close all
clc

load('N_disc_SS.mat', 'N_disc')

%% Options

% Remove nodes
remove = 1; % 0 - no nodes, 1 - remove section centred around (x,y)
x = 25; % Vertical location (k) of central removed node 
y = 5; % Horizontal location (j) of central removed node 
shape = 'horz'; % 'horz' - 1x3 section removed, 'vert' - 3x1 seciton removed

% Applied force 
applied = 1; %(N) - Ratio of applied force from standard loading conditions
E = 28E+09; %(Pa) - Young's Modulus of material

%% Parameters

T = 5; % Final time

% network size
K = 31; %Number of nodes (ROWS)
M = 11; %Number of nodes (COLS)

% signalling molecules
N_tot = 10000; % initial number of molecules 
N_ss = 400; % Expected number of molecules at steady state

Dt = 0.0005; % Time step interval

q_og = [0.05 0.05 0.05; 0.05 0.6 0.05; 0.05 0.05 0.05]; % Initial jump proportion matrix
q = q_og;

lambda1 = 100; % Generation of rate of molecules at each timestep   
lambda2 = 0.25; % Decay rate of molecules at each timestep 


% Force applied to lattice 
forces = applied* ones(M,1); % Uniform loading condiitons

%% Initialisation

t = 0; % Initial time 
iter = 0; % number of time steps taken
save_disc_every = 10; % Save solution state ever __ iterations
save_iter = 0; % counter: number of states saved
time_vec = [];

% Initalise flux
flux_left = 0;
flux_right = 0;
flux_up = 0;
flux_down = 0; 

% Initialise molecule counting vector 
total_part = zeros(T/Dt,1);

% Set reference psi value 
psi_bar = 0.5 * (1/E);

%% Time stepping 

% Establish data directory to store each assigned time step
mkdir('datadir')

save('datadir/state_disc_0.mat', 'N_disc', 'flux_left', 'flux_right', 'flux_up', 'flux_down') % save data for plotting later

N_disc_old = N_disc; 

N_disc_SS = N_disc; % keep a copy of SS for plotting contours. Will need existing N_disc saved in directory 

% Run FEM to calculate SED
if  remove == 1 
[SED] = simple_fem_remove(M,K, forces, applied, x, y, shape, psi_bar); % Remove nodes
else
[SED] = simple_fem(M, K, forces, applied, psi_bar); % No remove nodes
end 

% Time Stepping 
while t < T - 0.1*Dt % While the current time is less than the end time 
   
    % keep old values
    N_disc_old = N_disc;
    N_disc(2:K-1,2:M-1)=zeros(size(N_disc)-2);
  

        for j=2:M-1 % Loop over the lattice starting top left corner     
             for k=2:K-1

                if remove == 1 % Calculate node dependent q-matrix if nodes are removed 
                    N_disc(x,y) = 0; 
                    [q, lambda1_q, lambda2_q] = Newq_mult(j, k, x, y, q_og, lambda1, lambda2, shape);
                else q = q_og;
                    lambda1_q = lambda1;
                    lambda2_q = lambda2; 
                end 
    

             q = nofluxq(k, j, q, K, M); % Set q as no flux BC if along boundary (checked in script)


             % particles going out of j,k   
            for kk=-1:1
                for jj=-1:1
                    N_disc(k+kk, j+jj) = N_disc(k+kk, j+jj) + q(kk+2,jj+2) * N_disc_old(k,j);
                end
            end

            % Calculate molecule generation due to increased SED at each
            % element 

           psi = SED(k,j); % Extract SED for each node/element 
           N_disc(k,j) = N_disc(k,j) + Dt*(lambda1_q*(psi/ psi_bar) - lambda2_q*N_disc_old(k,j));
           
            end 
        end

    % Calculate flux
    flux_left = (N_disc(:,1) - N_disc_old(:,1))/Dt;
    flux_right = (N_disc(:,M) - N_disc_old(:,M))/Dt;
    flux_up = (N_disc(1,:) - N_disc_old(1,:))/Dt;
    flux_down = (N_disc(K,:) - N_disc_old(K,:))/Dt;
    
    % update time
    iter = iter + 1;
    t = t + Dt;


    % total particles
    total_part(iter) = sum(N_disc(:));

  

    % Save states to files state_1.mat, state_2.mat, etc.
    if rem(iter, save_disc_every) == 0
    
        save_iter = save_iter + 1;
        time_vec(save_iter) = t; 
        save(sprintf("datadir/state_disc_%d.mat", save_iter), 'N_disc', 'flux_left', 'flux_right', 'flux_up', 'flux_down')
    end

   
end

disp("Discrete model: done")

