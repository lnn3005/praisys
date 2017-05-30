%%%%%%%%%%% MRCPSP2_4.m
%%%%% as of 2016/10/04
%%%%% Multi-mode resource constrained project scheduling problem
%%%%% 2 modes, multiple(3) resource types
%%%%%
%%%%%

close all
clear
clc

%%%%%% read input
%%%%%%%% task matrix (task)
%%%%%%%% precedence matrix(pre) & successor matrix(suc)


%%% input file 

filename1 = 'input-task-42.txt';
filename2 = 'input-relation4.txt';

task = dlmread(filename1); 
rel = dlmread(filename2);   %%% precedence relations

%%%%% Ra: availability of renewable resources each period
%%%%%% how many resource types
m = 2; %%%% every task has 2 modes. 
mr = 3; %%%% # of resource type
Ra = [8 6 6];


%%% time horizon T = sum of durations of all tasks;

d1 = 2;
r1 = 3;
for i = 1:m;
    Dur(:,i)= task(:,d1+(i-1)*(mr+1));
    Res(:,mr*(i-1)+1:mr*i) = task(:,r1+((i-1)*(mr+1):(i-1)*(mr+1)+(mr-1)));
end



dur = [min(Dur,[],2), max(Dur,[],2)];  
tt = sum(dur);
T = max(tt); 




I = size(task,1); 



M = zeros(I,1)+m; M(1)=1; M(end)=1; %%% total number of modes

MM = zeros(length(rel),1);
for i = 1:length(MM)
    MM(i) = M(rel(i,1))* M(rel(i,2));
end



%%%%%%% precedence relationship
preMat = zeros(I,I);
for i = 1:length(rel)
    preMat(rel(i,1),rel(i,2)) = 1;
end

%%%%%% CPM5: activity on node (AON)
[elt,route,freetime,worktime] = CPM5(task,preMat,rel,dur,T);


%%%%%%%%%% construct matrices for intlinprog
%%%%% f
f = zeros(1,1+(T+1)*(I-2)*m+(T+1)); %%% 1+(T+1)*(I-2)*m+(T+1) = 1+(T+1)*((I-2)*m+1)
f((1+(T+1)*(I-2)*m+elt(end,3)+1):end) = linspace(elt(end,3),T,T-elt(end,3)+1);



%%%%%% intcon: index of integer location
intcon = linspace(1,1+(T+1)*((I-2)*m+1),1+(T+1)*((I-2)*m+1));


%%%%%%% inequality constraints
%%%%%%%%%%%%% A and b 
%%%%% A = [A1;A2]; A1-precedence, A2-resource constraint
%%%%% b = [b1;b2]; b1-precedence, b2-resource constraint


t = linspace(0,T,T+1);

%%%%%% A1 and b1
%%%%%% A1 shows the time LB and UB of tasks i and j, 
%%%%%% (i,j) belogns to precedence: i.e., task j can't starts until task i
%%%%%% finsihes. 
%%%%%% in A1: LBi=FEi,UBi = LFi, LBj=FEj,UBj = LFj




for i = 1:length(rel)
    tmpA1{i} = zeros(1,1+(T+1)*((I-2)*m+1));
    
    if MM(i) == 1
        tmpA1{i}(1,1) = 0;
        tmpA1{i}(1,(max(0,(rel(i,2)-2))*(T+1)*m+2)+(elt(rel(i,2),3):elt(rel(i,2),4))) = -t(1+elt(rel(i,2),3):1+elt(rel(i,2),4));

    end
    
    if MM(i) == m
        
        if rel(i,1) == 1 && rel(i,2) < task(end,1)
            tmpA1{i}(1,1) = 0;
            tmpA1{i}(1,(max(0,(rel(i,2)-2))*(T+1)*m+2)+(elt(rel(i,2),3):elt(rel(i,2),4))) = Dur(rel(i,2),1)-t(1+elt(rel(i,2),3):1+elt(rel(i,2),4));
            

            tmpA1{i}(1,(max(0,(rel(i,2)-2))*(T+1)*m+2+(T+1))+(elt(rel(i,2),3):elt(rel(i,2),4))) = Dur(rel(i,2),2)-t(1+elt(rel(i,2),3):1+elt(rel(i,2),4));
        
           
        else %if rel(i,1) > 1 && rel(i,2) == task(end,1)
            tmpA1{i}(1,(max(0,(rel(i,1)-2))*(T+1)*m+2)+(elt(rel(i,1),3):elt(rel(i,1),4))) = t(1+elt(rel(i,1),3):1+elt(rel(i,1),4));
            tmpA1{i}(1,(max(0,(rel(i,2)-2))*(T+1)*m+2)+(elt(rel(i,2),3):elt(rel(i,2),4))) = -t(1+elt(rel(i,2),3):1+elt(rel(i,2),4));
    
            tmpA1{i}(1,(max(0,(rel(i,1)-2))*(T+1)*m+2+(T+1))+(elt(rel(i,1),3):elt(rel(i,1),4))) = t(1+elt(rel(i,1),3):1+elt(rel(i,1),4));

        end
    
    end
    
    
    if MM(i) == m*m

     
        tmpA1{i}(1,(max(0,(rel(i,1)-2))*(T+1)*m+2)+(elt(rel(i,1),3):elt(rel(i,1),4))) = t(1+elt(rel(i,1),3):1+elt(rel(i,1),4));
        tmpA1{i}(1,(max(0,(rel(i,2)-2))*(T+1)*m+2)+(elt(rel(i,2),3):elt(rel(i,2),4))) = Dur(rel(i,2),1)-t(1+elt(rel(i,2),3):1+elt(rel(i,2),4));
         
        tmpA1{i}(1,(max(0,(rel(i,1)-2))*(T+1)*m+2+(T+1))+(elt(rel(i,1),3):elt(rel(i,1),4))) = t(1+elt(rel(i,1),3):1+elt(rel(i,1),4));
        tmpA1{i}(1,(max(0,(rel(i,2)-2))*(T+1)*m+2+(T+1))+(elt(rel(i,2),3):elt(rel(i,2),4))) = Dur(rel(i,2),2)-t(1+elt(rel(i,2),3):1+elt(rel(i,2),4));
        
   
    end
    
end


A1 = cell2mat(tmpA1');
b1 = zeros(length(rel),1);


%%%%%%%% Renewable resource constraints: A2 and b2
%%%%% For the time boudnary of renewable resource constraints
%%%%% LBi = max(t,EFi)
%%%%% UBi = min(t+djs-1,LFi), dis is the the duration of task i in mode s
%%%%% 

for it = 1:T+1
    for i = 1:I
        LBr(it,i) = max(t(it),elt(i,3));
        UBr_1(it,i) = min(t(it)+Dur(i,1)-1,elt(i,4));
        UBr_2(it,i) = min(t(it)+Dur(i,2)-1,elt(i,4));
    end
end
 

%%%%%%% resource availability in b2 vector
%%%% start with the 1st type of resource avialability
%%%% then the 2nd, and so forth
tmp = bsxfun(@plus,Ra,zeros(T+1,1));
b2 = reshape(tmp,numel(tmp),1);


%%%% find out locations where UNA2(it,i)>=LBA2(it,i), so that the summation
%%%% over time domain is doable. 
%%%% Relational Operation ge(A,B): Determine A is greater than or equal to B
%%%% then find row# and col# of nonzero elements in checkTimeBoundA2
checkTimeBoundr{1} = ge(UBr_1,LBr);
checkTimeBoundr{2} = ge(UBr_2,LBr);
[rnzo_1,cnzo_1] = find(checkTimeBoundr{1}~=0); 
[rnzo_2,cnzo_2] = find(checkTimeBoundr{2}~=0); 


   
%%%% for the 1st type of resource   
tmpA211 = zeros(1,1+(T+1)*((I-2)*m+1));   
   for i = 1:length(cnzo_1)
       if cnzo_1(i)>1 && cnzo_1(i)<task(end,1)
        
        tmpA211(rnzo_1(i),(max(0,(cnzo_1(i)-2))*(T+1)*m+2)+(LBr(rnzo_1(i),cnzo_1(i)):UBr_1(rnzo_1(i),cnzo_1(i)))) = Res(cnzo_1(i),1)+zeros(1,UBr_1(rnzo_1(i),cnzo_1(i))-LBr(rnzo_1(i),cnzo_1(i))+1);
            
       end
   end
   
tmpA212 = zeros(1,1+(T+1)*((I-2)*m+1));   
   for i = 1:length(cnzo_2)
       if cnzo_2(i)>1 && cnzo_2(i)<task(end,1) 
        
        tmpA212(rnzo_2(i),(max(0,(cnzo_2(i)-2))*(T+1)*m+2+(T+1))+(LBr(rnzo_2(i),cnzo_2(i)):UBr_2(rnzo_2(i),cnzo_2(i)))) = Res(cnzo_2(i),1+mr)+zeros(1,UBr_2(rnzo_2(i),cnzo_2(i))-LBr(rnzo_2(i),cnzo_2(i))+1);
            
       end
   end
   

%%%% for the 2nd type of resource     
tmpA221 = zeros(1,1+(T+1)*((I-2)*m+1));   
   for i = 1:length(cnzo_1)
       if cnzo_1(i)>1 && cnzo_1(i)<task(end,1)
        
        tmpA221(rnzo_1(i),(max(0,(cnzo_1(i)-2))*(T+1)*m+2)+(LBr(rnzo_1(i),cnzo_1(i)):UBr_1(rnzo_1(i),cnzo_1(i)))) = Res(cnzo_1(i),2)+zeros(1,UBr_1(rnzo_1(i),cnzo_1(i))-LBr(rnzo_1(i),cnzo_1(i))+1);
            
       end
   end
   
tmpA222 = zeros(1,1+(T+1)*((I-2)*m+1));   
   for i = 1:length(cnzo_2)
       if cnzo_2(i)>1 && cnzo_2(i)<task(end,1) 
        
        tmpA222(rnzo_2(i),(max(0,(cnzo_2(i)-2))*(T+1)*m+2+(T+1))+(LBr(rnzo_2(i),cnzo_2(i)):UBr_2(rnzo_2(i),cnzo_2(i)))) = Res(cnzo_2(i),2+mr)+zeros(1,UBr_2(rnzo_2(i),cnzo_2(i))-LBr(rnzo_2(i),cnzo_2(i))+1);
            
       end
   end

%%%% for the 3rd type of resource     
tmpA231 = zeros(1,1+(T+1)*((I-2)*m+1));   
   for i = 1:length(cnzo_1)
       if cnzo_1(i)>1 && cnzo_1(i)<task(end,1)
        
        tmpA231(rnzo_1(i),(max(0,(cnzo_1(i)-2))*(T+1)*m+2)+(LBr(rnzo_1(i),cnzo_1(i)):UBr_1(rnzo_1(i),cnzo_1(i)))) = Res(cnzo_1(i),3)+zeros(1,UBr_1(rnzo_1(i),cnzo_1(i))-LBr(rnzo_1(i),cnzo_1(i))+1);
            
       end
   end
   
tmpA232 = zeros(1,1+(T+1)*((I-2)*m+1));   
   for i = 1:length(cnzo_2)
       if cnzo_2(i)>1 && cnzo_2(i)<task(end,1) 
        
        tmpA232(rnzo_2(i),(max(0,(cnzo_2(i)-2))*(T+1)*m+2+(T+1))+(LBr(rnzo_2(i),cnzo_2(i)):UBr_2(rnzo_2(i),cnzo_2(i)))) = Res(cnzo_2(i),3+mr)+zeros(1,UBr_2(rnzo_2(i),cnzo_2(i))-LBr(rnzo_2(i),cnzo_2(i))+1);
            
       end
   end
   
   
A2 = [tmpA211+tmpA212;tmpA221+tmpA222;tmpA231+tmpA232];
   



% % A2 = [tmpA211+tmpA212];
% b2 = [];b2=tmp(:,1);
% 
% %%% for the 1st type of resource   
% tmpA21 = zeros(1,1+(T+1)*((I-2)*m+1));   
%    for i = 1:length(cnzo_1)
%        if cnzo_1(i)>1 && cnzo_1(i)<task(end,1)
%         
%         tmpA21(rnzo_1(i),(max(0,(cnzo_1(i)-2))*(T+1)*m+2)+(LBr(rnzo_1(i),cnzo_1(i)):UBr_1(rnzo_1(i),cnzo_1(i)))) = Res(cnzo_1(i),1)+zeros(1,UBr_1(rnzo_1(i),cnzo_1(i))-LBr(rnzo_1(i),cnzo_1(i))+1);
%             
%        end
%    end
%    
% tmpA22 = zeros(1,1+(T+1)*((I-2)*m+1));   
%    for i = 1:length(cnzo_2)
%        if cnzo_2(i)>1 && cnzo_2(i)<task(end,1) 
%         
%         tmpA22(rnzo_2(i),(max(0,(cnzo_2(i)-2))*(T+1)*m+2+(T+1))+(LBr(rnzo_2(i),cnzo_2(i)):UBr_2(rnzo_2(i),cnzo_2(i)))) = Res(cnzo_2(i),1+mr)+zeros(1,UBr_2(rnzo_2(i),cnzo_2(i))-LBr(rnzo_2(i),cnzo_2(i))+1);
%             
%        end
%    end
%        
% A2 = tmpA21+tmpA22;


A = [A1;A2];
b = [b1;b2];

%%%%% Equality constraints
%%%%% ensure that every task only executes once. 
%%%%%%%% Aeq and beq
%%%%%% time bounds of Aeq: LBi = EFTi, UBi = LFTi

beq = ones(I,1);



for i = 2:I-1 
   
   tmpAeq(i,(i-2)*(T+1)*m+1:(i-1)*(T+1)*m) = [zeros(1,elt(i,3)),ones(1,elt(i,4)-elt(i,3)+1),zeros(1,T-elt(i,4)+elt(i,3)),ones(1,elt(i,4)-elt(i,3)+1),zeros(1,T-elt(i,4))];
    
end

tmpAeq(I,((I-2)*(T+1)*m+1)+(elt(I,3):elt(I,4))) = ones(1,elt(I,4)-elt(I,3)+1);

Aeq = [zeros(I,1),tmpAeq];Aeq(1,1)=1;




%%%%%% lower bound and upper bound of X - decision variables
lb = zeros(1+(T+1)*((I-2)*m+1),1);
ub = ones(1+(T+1)*((I-2)*m+1),1);


tic
%%%%%%%%% intlinprog
[X,fval,exitflag] = intlinprog(f,intcon,A,b,Aeq,beq,lb,ub);



%%%% x1: the start dummay task is always done at t=0 in mode 1. 
%%%% x2: finish time at a chosen mode for task - 2:end-1 
%%%% Amode: find out which mode is chosen for every task in the schedule.

x1 = zeros(T+1,1);x1(1)=1;
x2 = reshape(X(2:end),T+1,(I-2)*m+1); 
x = [x1,x2];

mnzo = find(sum(x2)~=0);   
Amode = rem(mnzo,m); Amode(Amode==0)=m; Amode = [1,Amode]';
x2sum = x2(:,1:2:end) + [x2(:,2:2:end),zeros(T+1,1)]; x2sum = [x1,x2sum];


[r,c] = find(x2sum~=0); 
Ax = zeros(T+1,I);
for i = 1:I %%% actually length(r)=size(x2sum,2) = # of tasks = I
   
    Ax(r(i)-max(1,Dur(i,Amode(i)))+1:r(i),c(i)) = ones(max(1,Dur(i,Amode(i))),1);
    
end



for i = 1:mr
    rRes{i} = reshape([Res(2:end-1,i),Res(2:end-1,i+mr)]',(I-2)*m,1);
end
ArRes = cell2mat(rRes);

Ru = x(:,2:end-1)*ArRes;




figure, spy(Ax(2:fval+5,:)')
xlabel('time (period)')
ylabel('task #')

figure, bar(Ru(2:fval+5,:)) 
xlabel('time (period)')
ylabel('renewable resource used')
ylim([0 max(Ra)])

 






