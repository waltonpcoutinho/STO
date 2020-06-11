function [Cone,EndPlate1,EndPlate2] = Cone(X1,X2,R,n)
    % Calculating the height of the Cone
    height_cyl = norm(X2-X1);

    % Creating 2 circles in the XY plane
    theta = linspace(0,2*pi,n)';
    xa = R(1)*cos(theta);
    ya = R(1)*sin(theta);
    xb = R(2)*cos(theta);
    yb = R(2)*sin(theta);

    % Creating the points in the Z-Direction
    z1 = [0 height_cyl];

    % Creating (Extruding) the cylinder points in the X-Directions    
    X = [xa xb];
    Y = [ya yb];
    Z = repmat(z1,length(xa),1);

    Cone = surf(X,Y,Z);

    set(Cone,'XData',get(Cone,'XData')+X1(1));
    set(Cone,'YData',get(Cone,'YData')+X1(2));
    set(Cone,'ZData',get(Cone,'ZData')+X1(3));

    % Setting the color to the Cone and the end plates
    alpha = 0.5;
    set(Cone,'AmbientStrength',1,'FaceLighting','gouraud','FaceAlpha',alpha,'FaceColor','b');

    % Plot filled circle on the base of the cone
    hold on;
    fill(X1(1) + R(1)*cos(theta), X1(2) + R(1)*sin(theta),'b','FaceAlpha',alpha);
    
    % Plot border of the upper circle
    XX = X2(1) + R(2)*cos(theta);
    YY = X2(2) + R(2)*sin(theta);
    ZZ = Z(:,2);
    plot3(XX, YY, ZZ, 'Color', 'b');

    EndPlate1=[];
    EndPlate2=[];

    alpha2 = 0;
    set(Cone,'EdgeAlpha',alpha2);       
    
    axis equal;
    colormap jet;
    set(gca,'PlotBoxAspectRatio',[1 1 1]);
    grid on;
    hold on;
    
 end