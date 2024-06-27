 %% Flux along right boundary 

x_max = K-1; % Max x axis interval
x_min = 1; %Min x axis interval
Dxi = linspace(x_min, x_max, M);
Dyi = linspace(x_min, x_max, K);
s = 0;

% Select time steps of interest (can be any steps or length)
time_step = [0.005 0.01 0.05 0.1 0.3 0.5 1 2 3];

cmap = colormap(parula(length(time_step)));

for s=1:length(time_step)
    
    indx = time_step(s)/(Dt*save_disc_every); 
    load(sprintf("datadir/state_disc_%d.mat", indx), 'flux_left', 'flux_right', 'flux_up', 'flux_down') 
    plot(Dyi, flux_right, 'Color', [cmap(s,1), cmap(s,2), cmap(s,3)], 'LineWidth', 1.5)
    title('Flux at right boundary', 'FontSize', 12, 'FontName', 'times')
    xlabel('Position in osteon (osteocytes)', 'interpreter','latex','FontSize', 11, 'FontName', 'times')
    ylabel('Number of molecules', 'FontSize', 11, 'FontName', 'times')
    xticks([1 6 11 16 21 26 31] )
    xlim([0 31])
    ylim([0 10000]) %can be altered for scale
    hold on

end

legendCell = cellstr(num2str(time_step', 't = %.2f'));
leg = legend(legendCell , 'FontSize', 11, 'FontName', 'times');
leg.Location = 'eastoutside'; 

%print(gcf, sprintf(' *NAME FILE* - %f .png', applied),'-dpng','-r300')


%% Contour plot animation


f= figure;
cmap = colormap(parula(30));
m = linspace(0, 150, 30);

for s = -5:save_iter
    t = s*save_disc_every*Dt;

    if t <= 0 %plot initial steady state
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
    yticks([2 5 9 13 17 21 25 29])

     if  t == 5 %set given time to stop for saving image 

         print(gcf, sprintf(' *NAME FILE* %.3f - ratio %f .png', t, applied),'-dpng','-r300')

         break 
    end

    pause(0.0001)
    
end


%% Cross section plot 

figure 
Dxi = linspace(0, M, M);

% Select time steps of interest (can be any steps or length)
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
      
end


legendCell = cellstr(num2str(CS_times', 't = %.2f'));
leg = legend(legendCell , 'FontSize', 11, 'FontName', 'times');
leg.Location = 'eastoutside'; 

 print(gcf, sprintf(' *NAME FILE* - %f .png', applied),'-dpng','-r300')

