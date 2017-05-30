%%%%%%%%CPM5.m - unfinished
%%%% critical path method for activity on node (AON) network
%%% route = CPM without considering resource constraint
%%% worktime = final duration of CPM without considering resource constraint
%%% EFTj = earliest finish time of j
%%% LFTj = latest finish time of j


function [elt,route,freetime,worktime] = CPM5(task,preMat,relationship,dur,T)


% dur = task(:,2);
n = size(task,1);



for i = 1:n
   for j = 1:n
        if (preMat(i,j)==1)
            suc(i,j)=j;
            pre(j,i)=i;
        end
    end
end
        



for i = n-1:-1:1
    if (sum(suc(i,:)~=0,2)~=0)
        s = suc(i,find(suc(i,:)~=0));
    end
 
    ss(i,:) = [s, zeros(1,n-length(s))+s(end)];
    
end



est = zeros(n,1);
lst = zeros(n,1);
tmp1 = zeros(n,1);
for i = 2:n
    for j = 1:size(pre,2)
        if(pre(i,j)>0)
            tmp1(i,j) = dur(pre(i,j),1)+est(pre(i,j));
        end
    end
    est(i) = max(tmp1(i,:));
end


lst(n)=est(length(est));
tmp2 = zeros(n,1)+inf;


for i = n-1:-1:1
    for j = 1:size(ss,2)
       
        tmp2(i,j)=lst(ss(i,j))-dur(i,1);

    end

    lst(i) = min(tmp2(i,:));
end


eft = est + dur(:,1);
lft = zeros(size(task,1),1)+T;  lft(1) = 0;

% eft = est + dur(:,1);
% lft = lst + dur(:,2);


%%%%%% elt matrix = [est lst eft lft]
%%%%%%
elt = [est lst eft lft]; 

route=0;

for i=1:n
    if(est(i)==lst(i))
        route=[route i];
    end
end

for i=1:1:n
    freetime(i)=lst(i)-est(i);
end

route=route(2:length(route))';
worktime=est(n);


