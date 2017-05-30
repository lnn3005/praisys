### SETS

set V; # Set of tasks

set M{V}; # Set of modes

set R; # Set of resources

set P within {V,V}; # Precedence set

### PARAMETERS

param dummy;

param horizon_length;

set horizon:=1..horizon_length;

param e{V}; # Earliest start time

param l{V}; # Latest start time

param d{i in V, j in M[i]}; # Duration of task i in mode j

#param w{i in V, j in M[i], r in R}; # Resource required per period

param W{r in R,t in horizon}; # Available resources in period t

### DECISION VARIABLES

var x{i in V,j in M[i],horizon} binary;

var T{V};

var D{V};

var w{i in V, j in M[i], r in R};

### OBJECTIVE

minimize obj: T[dummy];

### CONSTRAINTS

st01 {i in V, j in M[i], r in R}: w[i,j,r] = 1;

st1 {i in V}: sum {j in M[i]} sum {t in e[i]..l[i]} x[i,j,t] = 1;

st2 {i in V}: T[i] =  sum {j in M[i]} sum {t in e[i]..l[i]} t * x[i,j,t];

st3 {i in V}: D[i] =  sum {j in M[i]} sum {t in e[i]..l[i]} d[i,j] * x[i,j,t];

st4 {(i,k) in P}: T[k] - T[i] >= D[i]; 

st5 {r in R, t in horizon}: sum {i in V} sum{j in M[i]} sum{theta in max(t,e[i])..min(l[i],t+d[i,j]-1)} w[i,j,r] * x[i,j,theta] <= W[r,t]; 

