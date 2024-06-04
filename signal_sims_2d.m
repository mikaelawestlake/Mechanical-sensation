
%% Create lattice pic 
h = 100/28;
vec = h/2 : h : 100-(h/2);
[X,Y] = meshgrid(vec,vec');
figure
plot(X, Y, 'k.', markersize=16)
axis image, axis([0, 100 ,0,100])
set(gca, 'clim', [0 100], 'YDir','reverse')
xlabel('j ($\mu m$)','Interpreter','latex')
ylabel('k ($\mu m$)','Interpreter','latex')

% use for exporting high res plot 
% print(gcf,'foo.png','-dpng','-r300');

%% 



if exist("datadir")==7
    rmdir("datadir", 's')
end

clear all
close all
clc

load('N_disc_SS.mat', 'N_disc')

%% Requirements


% Initialise time 
T = 5; % Final time

remove = 2;
x = 25;
y =5;
shape = 'horz'; %vertical stack or horizontal stack


% Impulse of molecules
imp = 0;
lam_imp = 0;
lam1_new = 10000;
lam2_new = 0.25;
lam_time = 100; 
lam_coord = [25 10];

%% Parameters

% network
K = 31; %Number of nodes (ROWS)
M = 11; %Number of nodes (COLS)

% signalling molecules
N_tot = 10000; % initial number of molecules 

N_ss = 400; % Expected number of molecules at steady state

% space and time discretisation
Dt = 0.0005; % Time step interval

q_og = [0.05 0.05 0.05; 0.05 0.6 0.05; 0.05 0.05 0.05]; 
q = q_og;

lambda1 = 100; 
lambda2 = 0.25;

applied = 1; %(N)
E = 28E+09; %(Pa)

% Force applied to lattice 
%forces = [zeros(1, M/2 + 0.5) applied*ones(1, M/2 -0.5)]; %step, force outside 
% forces = [applied*ones(1, M/2 -0.5) zeros(1, M/2 + 0.5)]; %step , force inside 
% forces = applied* linspace(-1,-2, M); % Distributed 
% forces = applied* linspace(-2,-1, M); % Distributed 
 forces = applied* ones(M,1); % Uniform
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

total_part = zeros(T/Dt,1);

psi_bar = 0.5 * (1/E)

%% Time stepping 

% Establish data directory to store each assigned time step
mkdir('datadir')

save('datadir/state_disc_0.mat', 'N_disc', 'flux_left', 'flux_right', 'flux_up', 'flux_down') % save data for plotting later

N_disc_old = N_disc;

N_disc_SS = N_disc; % keep a copy of SS for plotting contours 


if  remove == 1 || remove == 2 
[SED] = simple_fem_remove(M,K, forces, applied, x, y, shape, psi_bar);
else
 [SED] = simple_fem(M, K, forces, applied, psi_bar);
end 

while t < T - 0.1*Dt % While the current time is less than the end time 
   
    % keep old values
    N_disc_old = N_disc;
    N_disc(2:K-1,2:M-1)=zeros(size(N_disc)-2);
  

        for j=2:M-1 %Check bounds    
             for k=2:K-1
            if remove == 1
                N_disc(x,y) = 0; 
                [q, lambda1_q, lambda2_q] = Newq_mult(j, k, x, y, q_og, lambda1, lambda2, shape);
            elseif remove == 2 && strcmp(shape,'vert')  == 1      
                [q, lambda1_q, lambda2_q] = Newq_mult(j, k, x, y, q_og, lambda1, lambda2, shape);
            elseif remove == 2 && strcmp(shape,'horz')  == 1
                 [q, lambda1_q, lambda2_q] = Newq_mult(j, k, x, y, q_og, lambda1, lambda2, shape);
            else q = q_og;
                lambda1_q = lambda1;
                lambda2_q = lambda2; 
           
            end 
    

    % NO FLUX BC ON LOWER BOUNDARY 

             q = nofluxq(k, j, q, K, M); 

             % particles going out of j,k   

            for kk=-1:1
                for jj=-1:1
                    N_disc(k+kk, j+jj) = N_disc(k+kk, j+jj) + q(kk+2,jj+2) * N_disc_old(k,j);
                end
            end

            % lamb1 as funciton of SED
            psi = SED(k,j);
            
           N_disc(k,j) = N_disc(k,j) + Dt*(lambda1_q*(psi/ psi_bar) - lambda2_q*N_disc_old(k,j));
           
        end 
        end

    % fluxes
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

        if imp == 1
             if length(time_vec) == imp_time

                 N_disc(I_loc(1), I_loc(2)) = N_disc(I_loc(1), I_loc(2)) + impulse; 
    
             end
        end
 
    end

   
end

disp("Discrete model: done")




 %% Flux Map 


x_max = K-1; % Max x axis interval
x_min = 1; %Min x axis interval
Dxi = linspace(x_min, x_max, M);
Dyi = linspace(x_min, x_max, K);
s = 0;

time_step = [0.005 0.01 0.05 0.1 0.3 0.5 1 2 3];

cmap = colormap(parula(length(time_step)));

for s=1:length(time_step)
    
    indx = time_step(s)/(Dt*save_disc_every); 
    %indx = sprintf('%10.2f',indx)

    load(sprintf("datadir/state_disc_%d.mat", indx), 'flux_left', 'flux_right', 'flux_up', 'flux_down') 

    %subplot(2,3,2)
    plot(Dyi, flux_right, 'Color', [cmap(s,1), cmap(s,2), cmap(s,3)], 'LineWidth', 1.5)
    title('Flux at right boundary', 'FontSize', 12, 'FontName', 'times')
    xlabel('Position in osteon (osteocytes)', 'interpreter','latex','FontSize', 11, 'FontName', 'times')
    ylabel('Number of molecules', 'FontSize', 11, 'FontName', 'times')
    xticks([1 6 11 16 21 26 31] )
    xlim([0 31])
    ylim([0 10000])
    %xticklabels(linspace(0,10*K, length(xticks)))
    hold on

end

legendCell = cellstr(num2str(time_step', 't = %.2f'));
leg = legend(legendCell , 'FontSize', 11, 'FontName', 'times');
leg.Location = 'eastoutside'; 

 print(gcf, sprintf('horz out flux - %f .png', applied),'-dpng','-r300')


%% contour plot  


f= figure;

cmap = colormap(parula(30));

m = linspace(0, 150, 30);



for s = -5:save_iter
    t = s*save_disc_every*Dt;

    if t <= 0
    clf
    contourf(N_disc_SS(2:K-1, 1:M-1), m)
    title(sprintf( 'Time 0.000 '), 'FontSize',13, 'fontname','times');
else 
    load(sprintf("datadir/state_disc_%d.mat", s), 'N_disc') 
    clf
    contourf((N_disc(2:K-1, 1:M-1)), m)
    title(sprintf( 'Time %.3f', t), 'FontSize',13, 'fontname','times');
    end

    colormap(cmap);
    h = colorbar;
    set(get(h,'label'),'string','Number of molecules', 'Fontsize', 11, 'fontname','times') ;
    set(gca, 'clim', [0 100], 'YDir','reverse')
 
    
    ax.fontname = 'times';
    ay.fontname = 'times';
    axis equal
    xlabel('j', 'fontname','times', 'Fontsize', 11)
    ylabel('k', 'fontname','times', 'Fontsize', 11)
    xlim([2 M-1])
      xticks(2:10)
    % % xticklabels(linspace(0,5*M, leng(xticks)))
      yticks([2 5 9 13 17 21 25 29])
    % % yticklabels(linspace(0,35*K, length(xticks)))

    f.Position(3:4) = [300 500];

     if  t == 5

         print(gcf, sprintf('horz out - %.3f - ratio %f .png', t, applied),'-dpng','-r300')

         break 
    end

    pause(0.0001)
    
end



%% Cross section plot 

figure 

Dxi = linspace(0, M, M);

CS_times = [0.005 0.01 0.05 0.1 0.3 0.5 1 2 3];


cmap = colormap(parula(length(CS_times)));

for s = 1:length(CS_times)

    indx = CS_times(s)/(Dt*save_disc_every); 

    load(sprintf("datadir/state_disc_%d.mat", indx), 'N_disc') 
    to_plot = N_disc(x,:)';
    plot(Dxi(2:M-1), to_plot(2:M-1), 'Color', [cmap(s,1), cmap(s,2), cmap(s,3)],'LineWidth',1.5)
    ax.fontname = 'times';
    ay.fontname = 'times';
    
    xlabel('j', 'fontname','times', 'Fontsize', 11)
    ylabel('Number of molecules', 'fontname','times', 'Fontsize', 11)
    title(sprintf('Cross section of lattice over time'), 'FontSize',13, 'fontname','times');
    xlim([1 M-1])
    ylim([0 200])

    hold on 
         % need to re look over eq for analytical solution? cant find
         % source from pascal?? 
      
end

%  hold on 
% fplot(@(x)  0.3 * (lambda1/lambda2) * (1 - ( (cosh((sqrt(lambda2/100))* x)) / (cosh((sqrt(lambda2/100))*11)))), ...
%  'r --' , 'LineWidth',3   )

legendCell = cellstr(num2str(CS_times', 't = %.2f'));
leg = legend(legendCell , 'FontSize', 11, 'FontName', 'times');
leg.Location = 'eastoutside'; 

 print(gcf, sprintf('horz out CS - %f .png', applied),'-dpng','-r300')


%% 

fplot(@(x) 0.3 * (lambda1/lambda2) * (1 - ( (cosh((sqrt(lambda2/100))* x)) / (cosh((sqrt(lambda2/100))*11)))))
