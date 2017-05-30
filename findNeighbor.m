function [node,edge] = findNeighbor(nodes,links)

numNodes = size(nodes,1);
numEdges = size(links,1);
NodeSet = [1:numNodes];
EdgeSet = links;

% Choose root node. This can be interpreted as "depot" node
rootNode = 1; % -

% Choose root edge
rootEdgeIndex = 1;

% Determine the neighbor set for each node
for i=1:numNodes
    node(i).neighbor = [EdgeSet(EdgeSet(:,2) == i);EdgeSet(EdgeSet(:,1) == i,2)]';
end

% Determine the neighbor set for each edge
% All edges are referenced using INDICES, not pairs of nodes, e.g. edge 1 but not
% edge (1 2)
EdgeIndex = [1:numEdges]';
edge(rootEdgeIndex).neighbor = rootEdgeIndex; % Neighbor set of root edge is itself
for i=1:numEdges    
    if i == rootEdgeIndex
        continue;
    end
    
    edge(i).neighbor = cat(1,EdgeIndex(EdgeSet(:,1)==EdgeSet(i,1)),...
        EdgeIndex(EdgeSet(:,1)==EdgeSet(i,2)),....
        EdgeIndex(EdgeSet(:,2)==EdgeSet(i,1)),...
        EdgeIndex(EdgeSet(:,2)==EdgeSet(i,2)));
    edge(i).neighbor(edge(i).neighbor == i) = []; 
end

end