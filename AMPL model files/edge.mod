### SET ###

set N; #node set

set E within {N,N}; #edges set

param numPeriods; #number of time periods

set T = 1..numPeriods; #time horizon

set F{E} within E; #neighbor set of edge (i,j) in E

set G{N} within N; #neighbor set of node i in N

### PARAMETERS ###

param numCrewsNode; #number of available crews in each period

param numCrewsArc;

param workload{E}; #workload required on arcs

param functionality{E}; #initial functionality

param workloadNode{N};

### DECISION NARIABLES ###

var x{T,E} binary; #number of crews repairing arc e in period t 

var y{T,N} >=0, integer;

var u{T,E} binary; #accessibility of edge e at the START of period t

var v{T,E} binary; #functionality of edge e at the END of period t

var w{T,N} binary;

### OBJECTINE ###

maximize obj: sum {t in T, (i,j) in E} (i+j) * v[t,i,j];

### CONSTRAINTS ###

# st1 {t in T, (i,j) in E}: x[t,i,j] <= numCrewsArc*u[t,i,j];

st1 {t in T, (i,j) in E}: x[t,i,j] <= u[t,i,j];

st2 {t in T}: sum {(i,j) in E} x[t,i,j] <= numCrewsArc;

st3 {t in 2..numPeriods, (i,j) in E}: u[t,i,j] <= sum {(a,b) in F[i,j]} v[t-1,a,b];

st4 {t in T, (i,j) in E}: workload[i,j] * v[t,i,j] <= sum {r in 1..t} x[r,i,j];

st41 {(i,j) in E}: u[1,i,j] = functionality[i,j];

st5 {t in T, (i,j) in E}: v[t,i,j] <= u[t,i,j];

st6 {t in T, (i,j) in E}: x[t,i,j] >= 0;

# st7 {t in T, i in N}: y[t,i] <= numCrewsNode* sum{j in G[i]} v[t,i,j];

# st8 {t in T, i in N}: workloadNode[t,i]*w[t,i] <= sum {r in 1..t-1} y[r,i];

# st9 {t in T, i in N}: y[t,i] >= 0;
