
% the LINEAR PROGRAMMING APPROACH to DYNAMIC PROGRAMMING

n = 3; % states
m = 2; % actions

ggamma = 0.95; % discount factor

w = [1/3;1/3;1/3];

P{1} = [.8 .2 0
        .2 .6 .2
        .2 0 .8];  % transition matrix with action 1
    
P{2} = [.5 .5 0
        .3 .4 .3
        .7 0 .3];  % transition matrix with action 2
    
R1 = [2 -1 1
      -3 2 0
      4 0 -1];
  
R2 =  [1 0 2
       -2 1 1
        1 -2 1];

R{1} = sum(P{1}.*R1,2);  % rewards with action 1
R{2} = sum(P{2}.*R2,2);  % rewards with action 2



% cvx_solver gurobi
% cvx_solver mosek


% PRIMAL FORMULATION

cvx_begin
    variables v(n,1)  %value function
    dual variables llambda{m}
    minimize (w'*v)
    subject to
        for jj = 1:m
            llambda{jj}: v >= R{jj} + ggamma*P{jj}*v
        end
cvx_end

% DUAL FORMULATION

RR = [R{1}, R{2}];
cvx_begin
    variable x(n,m)
    dual variable vvv{n}
    maximize( trace(RR'*x) )
    subject to
     for s = 1:n
   vvv{s}:  sum(x(s,:)) - ggamma*[x(:,1)'*P{1}(:,s) + x(:,2)'*P{2}(:,s)] == w(s)
     end
    x >= 0
    
    sum(x(3,:)) >= 5
cvx_end

% OPTIMAL POLICY
pol = diag(1./sum(x,2))*x




% ROBUST MARKOV DECISION PROCESS
% assuming R is in a given ellipsoid
% R is minimized by an opponent 

RRbar = RR;
aalpha = 1.5;

cvx_begin
    variable x(n,m)
    dual variable vvv{n}
    maximize( trace(RRbar'*x) - aalpha*norm(x(:)))
    subject to
     for s = 1:n
   vvv{s}:  sum(x(s,:)) - ggamma*[x(:,1)'*P{1}(:,s) + x(:,2)'*P{2}(:,s)] == w(s)
     end
    x >= 0
cvx_end
pol_robust = diag(1./sum(x,2))*x


