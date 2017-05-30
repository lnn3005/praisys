function [nodes,links] = createLattice(numNodes)
% Create a grid-like network with numNodes nodes

side = floor(sqrt(numNodes));

numNodes = side^2;

% Columm vector
a = [1:side]';

% x-coordinate of nodes
x = repmat(a,side,1);

% y-coordinate of nodes
y = repmat(a,1,side)';

y = y(:);

% Nodes set
nodes = [0 0;x y];

% Number of links/edges
numLinks = 2*(side-1)*side+1;

% Initialize links/edges set
links = zeros(numLinks,2);
linkCounter = 2;
links(1,:) = [0 1];

% Create links/edges set
for i = 1:side
    for j = 1:side
        
        if j<side 
            links(linkCounter,:) = [side*(i-1)+j side*(i-1)+j+1];
            linkCounter = linkCounter+1;
        end
        
        if i<side 
            links(linkCounter,:) = [side*(i-1)+j side*i+j];
            linkCounter = linkCounter+1;
        end
        
    end
end


end
