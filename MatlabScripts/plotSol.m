function [] = plotSol(instance)

    labels = ["WORHP","IPOPT"];
    
    clc;
    fprintf("plotSol(%s)\n", instance);
    [inst,traj,minX,minY] = read(instance);
    fprintf('Plotting ... \n');
%     %plot waypoints
%     for i = 1: inst.noWaypts
%         plotCone(inst.waypts(i,:));
%     end
%     %plot landing sites
%     for i = 1: inst.noLandg
%         plotCircle(inst.landSites(i,:));
%     end
%     %plotTrajectory
%     lineHandle = [];
%     for k = 1: traj.nGliders        
%         gliderTrajectory = traj.gliders(k).trajectory(:,1:4);
%         nTimeSteps = traj.nTimeSteps;
%         routeSize = traj.routeSizes(k);
%         T = (routeSize - 1)*nTimeSteps;
%         lineHandle(k) = plotTrajectory(traj.type,gliderTrajectory,T,k);        
%     end
%     
%     legend(lineHandle, labels(1), labels(2));
%     
%     %set plot limits
%     ax = gca;
%     xlim(ax,[minX inf]);
%     ylim(ax,[minY inf]);
%     zlim(ax,[-10 inf]);
    
%     %plot state and control values
%     for k = 1: traj.nGliders
%         %compute time instants vector        
%         size = traj.routeSizes(k) - 1;
%         time = zeros(size*traj.nTimeSteps,1);
%         t0 = 0;
%         count = 1;
%         for i = 1: size
%             for j = 1: traj.nTimeSteps               
%                 time(count) = t0 + (j-1)*traj.gliders(k).steps(i);
%                 count = count + 1;
%             end            
%             t0 = time(count-1);
%         end       
%         gliderStatesAndControls = traj.gliders(k).trajectory(:,4:8);
%         plotVariables(gliderStatesAndControls, time, labels(k));
%     end
    
    %plot error analysis
    solver = ["IPOPT", "WORHP"];
    figure('Position', [0 0 1800 600]);
    tiledlayout(1,2);
    for k = 1: traj.nGliders
        nexttile;
        taylor_mean = mean(traj.gliders(k).maxNormTaylor);
        eps_mean = mean(traj.gliders(k).maxNormEps);
        aux = [taylor_mean; eps_mean]';
        bar(aux, 'stacked');
        hold on;
        
        yyaxis right
        relError_mean = mean(traj.gliders(k).relError);
        plot(relError_mean, 'r-', 'LineWidth', 3);
        
        ax = gca;
        ax.YAxis(2).Color = 'r';
        xlabel("Time steps");
        legend(ax, "Taylor's 1st term max .norm", ...
        "Epsilon max. norm", "Rel. error");
        title(solver(k));
    end
    
    %save figure to file
    figName = strcat(instance, '.jpg');
    saveas(gcf, figName);
  
end

