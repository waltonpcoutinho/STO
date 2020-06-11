function [] = plotVariables(variables,time)
    V = variables(:,1);
    pitch = variables(:,2);
    yaw = variables(:,3);
    CL = variables(:,4);
    roll = variables(:,5);  
        
    figure;
    hold on;
    grid on;
    plot(time,pitch);    
    plot(time,yaw);
    plot(time,roll);
    legend('Pitch','Yaw','Roll');
    xlabel('Time steps')
    
    figure;
    %steps = 1:1:length(V);    
    [hAx,~,~] = plotyy(time,V,time,CL);
    hAx(1).YColor = 'k';
    hAx(2).YColor = 'k';
    hold on;
    grid on;
    legend('V(t)','Cl(t)');
    xlabel('Time (s)');

end

