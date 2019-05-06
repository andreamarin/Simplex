%First Assignment.  
    %Andrea Marin Alarcon 158999
    %Andrea Perez Vega 154467

%The objective of the code is to solve (if it is possible) the following
%linear program:
%          maximise c^T x
%           subject to Ax = b, x >= 0, b >=0

%A = [2 1 0 1 0 0; 1 2 -2 0 1 0; 0 1 2 0 0 1];
%b = [10; 20; 5];
%c = [2; -1; 2; 0; 0; 0];

%A = [1 1 3 1 0; -2 0 2 0 1];
%c = [1; 1; 0; 0; 0];
%b = [5; -1];

%A = [4 1 1 0 0; 2 3 1 1 0; 1 2 3 0 1];
%b = [30; 60; 40];
%c = [3; 2; 1; 0; 0];

%A = [4 0 0 1 0 0; 6 1 0 0 1 0; 18 6 1 0 0 1];
%b = [1; 9; 81];
%c = [9; 3; 1; 0; 0; 0];

%A = [1 2 3/2 1 0 0;2/3 2/3 1 0 1 0; 1/2 1/3 1/2 0 0 1];
%b = [12000; 4600; 2400];
%c = [11; 16; 15; 0; 0; 0];

%A=[1 -1 1 0;-1 1 0 1];
%b=[1 2]';
%c=[1 0 0 0]';

A = [3 5; 4 1];
b = [78; 36];
c = [5;4];
[status, obasis, obfs, oval] = bothPhases(A,b,c)


%In this phase (2) we now want to find, it if exists, the optimal solution
%of the problem mentioned above, it solves also problems with complications
%(degeneracy).

% In the phase 2, we find the optimal solution of a LP problem, if it exists
% we do this by iterating over simplex tableaux using Bland's rule for the simplex method.
function[bound, obasis, obfs, oval] = phaseTwo(A, b, c, sbasis, sbfs)
% Input:
% A: mxn matrix with m <= n and rank of A is m.
% b: column vector with m rows.
% c: column vector with n rows.
% sbasis: vector of size m of indices of column vectors for a feasible basis for this
% problem from which to start the simplex method.
% sbfs: a vector of size n which is the basic feasible solution
% corresponding to the sbasis.

[m,n] = size(A);

A_B = zeros(m,m); % Nonsingular matrix that solves A_B*x_B=b; where B is the feasible basis and x_B the vector of the basic variables.
A_N = zeros(m, n-m);
c_B = zeros(m,1); 
c_N = zeros(n-m,1);

finished = false;

basic_var = sbasis;

while finished == false %This while is for the algorithm to repeat until it achieves the optimal solution.
    non_basic = setdiff(1:n, basic_var); %We obtain the nonbasic variables.

    A_B = A(:,basic_var);
    c_B = c(basic_var);
    A_N = A(:,non_basic);
    c_N = c(non_basic);

    Ab_inverse = inv(A_B);

    Q = -Ab_inverse * A_N; %Coefficients of the nonbasic variables in the functions of the basic variables.
    p = Ab_inverse * b; %Value of the basic variables when x_N=0.
    z_0 = transpose(c_B)* Ab_inverse * b; %Value of the objective function respect to the basic variables when x_N=0. 
    r = c_N - transpose(transpose(c_B) * Ab_inverse * A_N); %Coefficients of the nonbasic variables in the objective function (z).

    %We will use Bland's rule to avoid cycles.
    r_aux = find(r > 0); %We find the positive coefficients of the nonbasic variables in "z". 
    
    %If the coefficients of r are negative, then we are done.
    if isempty(r_aux)
        finished = true;
        bound = 1; %The problem is bounded, and therefore, it exists an optimal solution.
        continue
    end
    
    % We obtain the indices of the variables that have positive
    % coefficients in r.
    indices = zeros(1, length(r_aux));
    for i = 1:1:length(r_aux)
        index = r_aux(i);
        indices(i) = non_basic(index);
    end
    enters = min(indices);
    enters_pos = find(non_basic == enters);
    
    % We develop the condition that the leaving variable has to satisfy.
    q_aux = find(Q(:,enters_pos) < 0);
    
    %If there are no negative inputs in Q, the problem is not bounded,
    %therefore you finish.
    if isempty(q_aux)
        finished = true;
        bound = 0; %The problem is unbounded and it doesn't have an optimal solution.
        continue
    end
    
    constraints = zeros(1, length(q_aux));
    for i = 1:1:length(q_aux)
        index = q_aux(i);
        constraints(i) = -p(index)/Q(index,enters_pos);
    end

    %We find the minimum constraint for the leaving variable.
    min_constraint = min(constraints);
    pos = find(constraints == min_constraint);
    
    if length(pos) > 1
        indices = zeros(length(pos),1);
        for i = 1:1:length(pos)
            indices(i) = basic_var(q_aux(pos(i)));
        end
        leaves = min(indices);
    else
        leaves = basic_var(q_aux(pos(1)));
    end
    
    index = find(basic_var == leaves);
    basic_var(index) = enters;
    
end

oval = z_0; %the objective value of this optimal basic feasible solution (if the problem
% is bounded)
obfs = zeros(1,n); %vector of size n which is the optimal basic feasible solution 
% corresponding to the obasis.
obasis = basic_var; %vector of size m of indices of column vectors which gives an optimal fea-
% sible basis for the problem if the problem is bounded.

for i = 1:1:m
    index = basic_var(i);
    obfs(index) = p(i);
end

end

