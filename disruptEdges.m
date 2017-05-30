function [init_acc,workload] = disruptEdges(nodes,links)
% Create disruption on edges of a lattice network.
% Functionality is a binary vector representing initial accessibility,
% which is vector u_1 in model.

side = floor(sqrt(size(nodes,1)));
numLinks = 2*(side-1)*side+1;

% Destroy some links/edges
% functionality = randi(100,[numLinks 1]) > percentage;
init_acc = zeros([numLinks 1]);
init_acc(1) = 1;
init_acc(2) = 1;
init_acc(3) = 1;

workload = randi(5,[numLinks 1]) - 1;

% Workload of root edge must be 0
workload(1) = 0;

% Workload of edge 2 and 3 must be non-zero. Edge 2 and 3 are in the
% neighbor set of root edge.
workload(2) = 1;
workload(3) = 1;

% edge = findNeighbor(nodes,links);
% 
% for i = 2:numLinks
%     
%     if functionality(i) == 1
%         continue;
%     end
%     
%     if sum(functionality(edge(i).neighbor)) == 0
%         workload(i) = workload(i)*(randi(100) > 20);
%     end
%     
% end

end


