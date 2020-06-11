function [] = plotCone(waypoint)
    %waypoint parameters
    basePosition = waypoint(1:3);
    topPosition = basePosition;
    topPosition(3) = topPosition(3) + waypoint(6);

    baseRadius = waypoint(4);
    h = baseRadius/tan(45*pi/180);
    H = topPosition(3);
    topRadius = ((H + h)*baseRadius)/h;

    Radii = [baseRadius topRadius];

    %fixed parameters
    n=100;
    cyl_color='b';
    closed=0;
    lines=0;
    alpha=0.2;

    %plot one cone
    %[Coneh,EndPlate1,EndPlate2] = Cone(X1,X2,R,n,cyl_color,closed,lines,alpha);
    [~,~,~] = Cone(basePosition,topPosition,Radii,n);
end

