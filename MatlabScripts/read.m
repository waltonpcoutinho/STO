function [instData,trajData,minX,minY] = read(instance)   
    %read instance
    addpath('../../Python/InstanceGenerator/S/');
    addpath('../../Python/InstanceGenerator/M/');
    addpath('../../Python/InstanceGenerator/L/');
    addpath('../../Python/InstanceGenerator/Real_instances/');
    
    instFile = strcat(instance,'.dat');
    [instData] = readInstData(instFile);
    fprintf('Instance file read!\n');
        
    %read optimal trajectory
    solFile = strcat('../Results/Solutions/solution_',instance,'.out');
    [trajData, minX, minY] = readTrajData(solFile);
    fprintf('Solution file read!\n');

end

function [instData] = readInstData(instFile)
    %open file
    fileID = fopen(instFile);
    %read depot coordinates
    depotCoord = str2num(fgetl(fileID));
    %read vertices coordinates
    noVertices = str2num(fgetl(fileID));    
    for i = 1: noVertices
        waypoints(i,:) = str2num(fgetl(fileID));
    end
    %read landing points coordinates    
    noLanding = str2num(fgetl(fileID));    
    for i = 1: noLanding
        landing(i,:) = str2num(fgetl(fileID));
    end
    %export data to struct and close file        
    instData = struct('depot',depotCoord,'noWaypts',noVertices,'waypts',waypoints,'noLandg',noLanding,'landSites',landing);
	fclose(fileID);
end

function [trajData, minX, minY] = readTrajData(name)
    %open file
    fileID = fopen(name);
    %read instance name
    instance = fgetl(fileID);
    %find type of instance
    if (instance(6) == 'S')
        type = 'small';
    elseif (instance(6) == 'M')
        type = 'medium';
    elseif (instance(6) == 'L')
        type = 'large';
    end
    %read number of time steps
    nTimeSteps = str2num(fgetl(fileID));
    %read number of gliders
    nGliders = str2num(fgetl(fileID));
    %get the second element which is the # of used gliders
    nGliders = nGliders(2);  
    %read vector of route sizes
    routeSizes = str2num(fgetl(fileID));
    %read routes
    for i = 1: nGliders
        gliders(i).route = str2num(fgetl(fileID));
    end
    %read flight times of each route
    for i = 1: nGliders
        gliders(i).flightTimes = str2num(fgetl(fileID));
    end
    %read step sizes of each route
    for i = 1: nGliders
        gliders(i).steps = str2num(fgetl(fileID));
    end
    %read errors of each route
    for i = 1: nGliders
        gliders(i).errors = str2num(fgetl(fileID));
    end
    
    %read glider's positions
    for i = 1: nGliders
        for j = 1: (routeSizes(i) - 1)*nTimeSteps            
            line = fgetl(fileID);
            if(ischar(line))
                trajectory(j,1:8) = str2num(line);
            else
                fprintf('Error reading solution file!')
                return;
            end
        end
        gliders(i).trajectory = trajectory;
        trajectory = [];
    end
    
    minX = 9999999999;
    minY = 9999999999;
    
    for i = 1: nGliders
        minXAux = min(gliders(i).trajectory(1:end,1));
        if minXAux < minX
            minX = minXAux;
        end
        
        minYAux = min(gliders(i).trajectory(1:end,2));
        if minYAux < minY
            minY = minYAux;
        end
    end
    
    trajData = struct('instName',instance,'type',type,'nTimeSteps',nTimeSteps,'nGliders',nGliders,...
                      'routeSizes',routeSizes,'gliders',gliders);
	fclose(fileID);
end

