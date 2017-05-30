function [nodes,links] = createTransportationNetwork(numNodes) 

distanceWeight = 6;

nodes = randi(numNodes*distanceWeight,[numNodes 2]);

% Remove duplicated nodes and sort positions.
nodes = unique(nodes,'rows');

nodes = nodes(1:distanceWeight:size(nodes,1),:);

scatter(nodes(:,1),nodes(:,2));




end
