/* TCP, Traveling Couples Problem */

/* number of nodes */
param n, integer >= 3;

/* distance in meters */
param d, integer >= 0;

/* Set of nodes */
set V := 1..n;

/* Set of edges */
set E, within V cross V;

/* Distance */
param c{(i,j) in E};

/* x[i,j] = 1 means the edge is in the path */
var x{(i,j) in E}, binary;

/* Objective function */
minimize total: sum{(i,j) in E} c[i,j] * x[i,j];

/* Leaves once or not at all */
s.t. leave{i in V}: sum{(i,j) in E} x[i,j], <= 1;

/* Enters once or not at all */
s.t. enter{j in V}: sum{(i,j) in E} x[i,j], <= 1;

/* enter = leave */
s.t. samedef{i in V}: sum{(i,j) in E} x[i,j] - sum{(j,i) in E} x[j,i], = 0;

/* Vicinity */
s.t. vicinity{k in V}:
    sum{(i,j) in E} (if (c[i,k] <= d or c[j,k] <= d or i = k or j = k)
    then x[i,j] else 0), >= 2;

/* subtour stuff: */
var y{(i,j) in E}, >= 0;

/* Upper bound for y[i,j]. Mainly used so that (i,j) not in path = 0 */
s.t. cap{(i,j) in E}: y[i,j] <= (n-1) * x[i,j];

s.t. subtour{i in V}:
    sum{(j,i) in E} y[j,i]
    + (if i = 1 then sum{(j,k) in E} x[j,k])
    =
    sum{(i,j) in E} y[i,j]
    + 1 * (sum{(i,j) in E} x[i,j]);

solve;

printf "Optimal tour has length %d\n",
    sum{(i,j) in E} c[i,j] * x[i,j];
printf("From node   To node   Distance\n");
printf{(i,j)in E: x[i,j]} "     %3d     %d    %8g\n",
    i, j, c[i,j];
/*
data;

param n := 9;
param d := 35;

param : E : c :=
    1 9 173
    1 2 164
    1 3 184
    1 4 193
    1 5 184
    1 6 165
    1 7 154
    1 8 174
    1 1 10000000
    9 2  20
    9 3  20
    9 4  20
    9 5  20
    9 6  20
    9 7  20
    9 8 173
    9 1 173
    9 9 1000000
    2 9  20
    2 3  20
    2 4  34
    2 5  40
    2 6  34
    2 7  20
    2 8 184
    2 1 164
    2 2 1000000
    3 9  20
    3 2  20
    3 4  20
    3 5  40
    3 6  34
    3 7  20
    3 8 193
    3 1 184
    3 3 1000000
    4 9  20
    4 2  34
    4 3  20
    4 5  20
    4 6  34
    4 7  39
    4 8 144
    4 1 193
    4 4 1000000
    5 9  20
    5 2  40
    5 3  34
    5 4  20
    5 6  20
    5 7  34
    5 8 164
    5 1 184
    5 5 1000000
    6 9  20
    6 2  34
    6 3  39
    6 4  34
    6 5  20
    6 7  20
    6 8 154
    6 1 165
    6 6 1000000
    7 9  20
    7 2  20
    7 3  34
    7 4  39
    7 5  34
    7 6  20
    7 8 165
    7 1 154
    7 7 1000000
    8 9 173
    8 2 184
    8 3 193
    8 4 184
    8 5 164
    8 6 154
    8 7 165
    8 1 174
    8 8 1000000
;
*/
end;
