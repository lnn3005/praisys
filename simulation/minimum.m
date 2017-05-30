function result = minimum(nums,m)
% Return indices of m minimum numbers in vector nums
result = zeros(1,m);

if size(nums,1) > 1
    nums = nums';
end

if size(nums,1) > 1 
    display('Input must be a vector');
    return;
end

copy = nums;
counter = 1;
while counter <= m
    min = Inf;    
    for i=1:size(nums,2)
        if min>nums(i)
            min = nums(i);
            result(counter) = i;
        end
    end
    nums(result(counter)) = Inf;
    counter = counter+1;
end

% Restore content of nums
nums = copy;

clear copy;
end