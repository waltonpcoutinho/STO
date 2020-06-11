function [] = plotSol(instance)
    clc;
    fprintf("plotSol('grtopL_103_2')\n");
    [inst,traj,minX,minY] = read(instance);
    fprintf('Plotting ... \n');
    %plot waypoints
    for i = 1: inst.noWaypts
        plotCone(inst.waypts(i,:));
    end
    %plot landing sites
    for i = 1: inst.noLandg
        plotCircle(inst.landSites(i,:));
    end
    %plotTrajectory
    lineHandle = [];
    for k = 1: traj.nGliders        
        gliderTrajectory = traj.gliders(k).trajectory(:,1:4);
        nTimeSteps = traj.nTimeSteps;
        routeSize = traj.routeSizes(k);
        T = (routeSize - 1)*nTimeSteps;
        lineHandle(k) = plotTrajectory(traj.type,gliderTrajectory,T,k);        
    end
    
    legend(lineHandle, 'STO','STO-NLP');
    
    %set plot limits
    ax = gca;
    xlim(ax,[minX inf]);
    ylim(ax,[minY inf]);
    zlim(ax,[-10 inf]);

    %save figure to file
    %figName = strcat(instance, '.fig');
    %savefig(figName);
  
end

