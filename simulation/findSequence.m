function [sequence,time] = findSequence(x)
% Given soluion matrix x, this function return the sequence of fixing edges
% and the total time spent fixing them.

n = size(x,2);
sequence = zeros(1,n);

% Total time spent on fixing each edge
time = sum(x);

counter = 1;
for i = 1:size(x,1)
    for j = 1:n
        if x(i,j) == 1 && sequence(j) == 0
            sequence(j) = counter;
            counter = counter + 1;
        end
    end
end


