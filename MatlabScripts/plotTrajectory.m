function [lineHandle] = plotTrajectory(instType,singleTraj,T,index)

    % add path to 3D airplane model
    addpath('c130_v1/');
    
    X = singleTraj(:,1);
    Y = singleTraj(:,2);
    Z = singleTraj(:,3);
    V = singleTraj(:,4);
    
    %plot central trajectory
    lineWidth = 1;
    
    %glider body colour
    color = 'black';   
    %glider tailwing colour
    tailwing = {'blue'; 'red'; 'green'; 'magenta'; 'cyan'; [1.0 0.4 0.0]; [0.8 0.0 1.0]; [0.0 1.0 0.6]; [0.6 0.0 0.0]; [0.0 0.2 0.0]; [0.0 0.0 0.2]; [0.8 0.0 0.0]; [0.6 0.0 1.0]; [0.8 0.0 0.6]};
    
    %approximate central trajectory by CUBIC splines
    %note that the input of cscvn is TRANSPOSED
    smoothTrajectory = [X Y Z]';
    func = cscvn(smoothTrajectory(:,[1:end]));   
    points = fnplt(func, 2)';    
    lineHandle = plot3(points(:,1),points(:,2),points(:,3),'LineWidth',lineWidth,'Color',tailwing{index});
    lineHandle = lineHandle(1);
    hold on;
    
    % approximate angles from trajectory
    [mu,gamma,phi] = xyz2rpy(points(:,1),points(:,2),points(:,3));

    %define scale of model glider
    if (strcmp(instType,'small'))
        scale = 2;
    elseif (strcmp(instType,'medium'))
        scale = 5;
    elseif (strcmp(instType,'large'))
        scale = 10;
    end
       
    %plot UAVs
    for t = 1:20:length(points)
        handle = c130(points(t,1),points(t,2),points(t,3),'color',color,'tailwing',tailwing{index},...
               'lines','none','pitch',gamma(t),'yaw',phi(t),'roll',mu(t),...
               'scale',scale);
        hold on;
    end

    xlabel('$x$','Interpreter','LaTex');
    ylabel('$y$','Interpreter','LaTex');
    zlabel('$h$','Interpreter','LaTex');
        
    grid on;
    hold on;
end

