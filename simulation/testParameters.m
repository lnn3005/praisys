%% Run simulation for different parameters

% Minimum and maximum length of sides in network
minSide = 3;
maxSide = 50;

% Minimum and maximum number of crews
minCrews = 1;
maxCrews = 20;

% Rows and columns of result matrix
rows = maxSide - minSide + 1;
cols = maxCrews - minCrews + 1;

% Initialize result matrices
time = zeros(rows,cols);
serviceLevel = zeros(rows,cols);
average_func = zeros(rows,cols);

% Maximum iterations
maxIteration = 5;

for i = minSide:maxSide 
    sizes = i^2;
    numEdges = 2*i*(i-1);
    
    for j = minCrews:maxCrews
        crews = j;
        
        a = zeros(maxIteration,1);
        b = zeros(maxIteration,1);
        c = zeros(maxIteration,1);
        for k = 1:maxIteration  
            [a(k),b(k),c(k)] = runSimulation(sizes,crews);
        end
        
        time(i,j) = sum(a)/maxIteration;
        serviceLevel(i,j) = sum(b)/maxIteration - 48;
        average_func(i,j) = sum(c)/maxIteration;
    end
    
end

X = minSide:maxSide;
Y = minCrews:maxCrews;

figure;
surf(time(X,Y));
title('Average time required to recover network');
xlabel('Number of crews');
ylabel('Size of network');

figure;
surf(average_func(X,Y));
title('Average functionality of network');
xlabel('Number of crews');
ylabel('Size of network');

figure;
surf(serviceLevel(X,Y));
title('Service level');
xlabel('Number of crews');
ylabel('Size of network');

