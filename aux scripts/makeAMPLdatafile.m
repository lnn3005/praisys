%% Get input

% Hardcoded for now.

numPeriods = 7;

numCrews = 7;

% Network input
run('input_transportation_example1.m');

numNodes = size(ST.nodes,1);
numEdges = size(ST.links,1)/2;
NodeSet = [1:numNodes];
EdgeSet = ST.links(:,1:2);
EdgeSet1 = EdgeSet(1:2:size(EdgeSet,1),:); % No-duplication edge set


%% Neighbor works

% Choose root node. This can be interpreted as "depot" node
rootNode = 1; % -

% Choose root edge
rootEdgeIndex = 7;

% Determine the neighbor set for each node
for i=1:numNodes
    node(i).neighbor = EdgeSet(EdgeSet(:,2) == i)';
end

% Determine the neighbor set for each edge
% All edges are referenced using indices, not nodes, e.g. edge 1 but not
% edge (1 2)
EdgeIndex = [1:numEdges]';
edge(rootEdgeIndex).neighbor = rootEdgeIndex; % Neighbor set of root edge is itself
for i=1:numEdges    
    if i == rootEdgeIndex
        continue;
    end
    
    edge(i).neighbor = cat(1,EdgeIndex(EdgeSet1(:,1)==EdgeSet1(i,1)),...
        EdgeIndex(EdgeSet1(:,1)==EdgeSet1(i,2)),....
        EdgeIndex(EdgeSet1(:,2)==EdgeSet1(i,1)),...
        EdgeIndex(EdgeSet1(:,2)==EdgeSet1(i,2)));
    edge(i).neighbor(edge(i).neighbor == i) = []; 
end

%% Additional parameters

% Workload
workload = [0 0 4 0 5 3 0 7 6 8 4 5];

% Functionality
functionality = [1 1 0 0 0 0 1 0 0 0 0 0];

workloadString = [EdgeSet1';workload];
workloadString = workloadString(1:numEdges*3);

funcString = [EdgeSet1';functionality];
funcString = funcString(1:numEdges*3);

%% Start AMPL data file 

% Create a subfolder called 'AMPL data files'. A new data file is created
% inside this folder each time the code is run.

dir = 'AMPL data files';
filename = strcat(datestr(now,'ddmmmyy_HHMM'),'.dat');

if ~exist(dir,'dir')  
    mkdir(dir);
end

%% Write to .dat file

% Open data file to be written
fileID = fopen(fullfile(dir,filename),'w');

fprintf(fileID,'data; \n\n');

% Print node set
writeAMPL(fileID,1,'N',NodeSet);

% Print edge set
writeAMPL(fileID,11,'E',EdgeSet1);

% Print neighbor node sets
%{
for i=1:numNodes
    nameString = ['F[',num2str(i),']'];
    writeAMPL(fileID,1,nameString,node(i).neighbor);
end
%}

% Print neighbor edge sets
for i=1:numEdges
    nameString = ['F[',num2str(EdgeSet1(i,1)),',',num2str(EdgeSet1(i,2)),']'];
    writeAMPL(fileID,11,nameString,EdgeSet1(edge(i).neighbor,:));
end

% Print parameter
writeAMPL(fileID,2,'numPeriods',numPeriods);
writeAMPL(fileID,2,'numCrews',numCrews);

writeAMPL(fileID,21,'workload',workloadString);
writeAMPL(fileID,21,'functionality',funcString);


fclose(fileID);