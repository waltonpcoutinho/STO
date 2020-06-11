function [] = plotCircle(landing)
    %x and y are the coordinates of the center of the circle
    %r is the radius of the circle
    %0.01 is the angle step, bigger values will draw the circle faster but
    %you might notice imperfections (not very smooth)
    
    %fixed parameters
    n=50;
    red_color='r';
    alpha=0.2;
    
    %landing point data    
    x = landing(1);
    y = landing(2);
    z = landing(3);
    r = landing(4);
    
    phi = linspace(pi/2,-pi/2,n);
    theta = linspace(pi/2,-pi/2,n);
    [phi,theta] = meshgrid(phi,theta);

    xp = r*sin(phi).*cos(theta);
    yp = r*sin(phi).*sin(theta);
    zp = r*cos(phi); 
        
    %plot landing site
    surf(x + xp, y + yp, z + zp,'AmbientStrength',1,'FaceColor',red_color,...
    'FaceLighting','gouraud','FaceAlpha',alpha,'EdgeAlpha',0);%, 'MarkerEdgeColor','r');
    grid on;
    hold on;
end
