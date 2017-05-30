%% INPUT

[nodes,links] = createLattice(9);
% nodes is the co-ordinate for each node. links is the edges indexed by node.

[node,edge] = findNeighbor(nodes,links);
% Both node and edge are structs

[init_acc,workload] = disruptEdges(nodes,links);


numNodes = size(nodes,1);
numEdges = size(links,1);

NodeSet = [0:numNodes-1];
EdgeSet = links;

% workloadString = [EdgeSet';workload'];
% workloadString = workloadString(1:numEdges*3);

% funcString = [EdgeSet';functionality'];
% funcString = funcString(1:numEdges*3);

numPeriods = 9;

numCrewsEdges = 5;

% workload = [0 1 1 4 4 3 0 1 1 3 0 3 0]';

%% Deterministic simulation

% Vector of remaining workload
% remain = workload';
remain = (workload > 0)';

% Initial accessibility matrix u
u = zeros(numPeriods,numEdges);
u(1,:) = init_acc;

% Initial functionality matrix v
v = zeros(numPeriods,numEdges);
v(:,1) = 1; % Root edge must be functional

% Action matrix x
x = zeros(numPeriods,numEdges);

t = 1;
finish_in_time = 'No';
time = numPeriods;
% Main loop - stop when all edges are functional
while sum(v(t,:))<numEdges && t<numPeriods
    
    availCrews = numCrewsEdges - sum(x(t,:));
    
    % Update u - accessibility
    if t>1
        % Loop through all edges
        for i = 2:numEdges
            % Update u if the edge is inaccessible and the sum of
            % functionality of neighbor egdes is positie
            if u(t,i) == 0 && sum(v(t,edge(i).neighbor))>0 
                % Update u for all subsequent periods
                u(t:numPeriods,i) = 1;
            end            
        end   
    end
    
    % Compute vector of possible actions - a set of candidate edges
    actionSpace = u(t,:).*(v(t,:)==0).*(x(t,:)==0).*(remain>0);
    
    % Maximum number of crews that can be allocated in this period
    maxCrew = min(availCrews,sum(actionSpace));
    
    % newWorkload vector contains only the workload of candidate edges.
    % Workload on other edges is set to Inf
    newWorkload = workload'.*actionSpace;
    newWorkload(newWorkload==0) = Inf;

    % Assign crews to edges
    x(t,minimum(newWorkload,maxCrew)) = 1;
    
    % Decide if the chosen edges can be completely recover in this period
    for i=2:numEdges
        if x(t,i) == 1
            remain(i) = rand() > 1/workload(i);
            if remain(i) > 0
                x(t+1,i) = 1;
                % fprintf('Edge %d in period %d is not fixed \n',i,t);
            end
        end
    end 
    
    % Update v - functionality
    for i = 2:numEdges
        % Update v(i) if edge i is accessible and remaining workload on i
        % is less than or equal to 0.
        if v(t,i) == 0 && u(t,i) == 1 && remain(i) <= 0
            v(t:numPeriods,i) = 1;            
            % Free up a crew
        end
    end
    
    if sum(v(t,:)) == numEdges
        finish_in_time = 'Yes';
        time = t;
    end

    % Increment t
    t = t + 1;
end

fprintf('Finish in time? %s \nTime = %d \n',finish_in_time, time);

