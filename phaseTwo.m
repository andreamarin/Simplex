% First Assignment.  
%   Andrea Marin Alarcon 158999
%   Andrea Perez Vega 154467
%    
% The objective of the code is to solve (if it is possible) the following
% linear program:
%          maximise c^T x
%           subject to Ax = b, x >= 0, b >=0


% In Phase Two, we find the optimal solution of a LP problem, if it
% exists.
% We do this by iterating over simplex tableaux using Bland's rule for the simplex method.
function[bound, obasis, obfs, oval] = phaseTwo(A, b, c, sbasis, sbfs)
% INPUT:
% A: mxn matrix with m <= n and rank of A is m.
% b: column vector with m rows.
% c: column vector with n rows.
% sbasis: vector of size m of indices of column vectors for a feasible basis for this
% problem from which to start the simplex method.
% sbfs: a vector of size n which is the basic feasible solution
% corresponding to the sbasis.

% OUTPUT:
% bound = 0 if the problem is unbounded, 1 otherwhise
% obasis = basis for the optimal solution
% obfs = optimal basic feasible solution
% oval = optimal value for the objective function

[~,n] = size(A);
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
    indices = non_basic(r_aux);
    %indices = zeros(1, length(r_aux));
    %for i = 1:1:length(r_aux)
    %    index = r_aux(i);
    %    indices(i) = non_basic(index);
    %end
    
    enters = min(indices);
    enters_pos = find(non_basic == enters);
    
    % We develop the condition that the leaving variable has to satisfy.
    q_aux = find(Q(:,enters_pos) < 0);
    
    % If there are no negative inputs in Q, the problem is not bounded,
    % therefore you finish.
    if isempty(q_aux)
        finished = true;
        bound = 0; %The problem is unbounded and it doesn't have an optimal solution.
        continue
    end
    
    % find the entering variable's constraints
    constraints = -p(q_aux)/Q(q_aux,enters_pos);
    
    %constraints = zeros(1, length(q_aux));
    %for i = 1:1:length(q_aux)
    %    index = q_aux(i);
    %    constraints(i) = -p(index)/Q(index,enters_pos);
    %end

    %We find the minimum constraint for the leaving variable.
    min_constraint = min(constraints);
    min_pos = q_aux(constraints == min_constraint);
    
    % if there is more than one candidate for the leaving variable we
    % choose the one with the smallest index
    if length(min_pos) > 1
        indices = basic_var(min_pos);
        
        %indices = zeros(length(pos),1);
        %for i = 1:1:length(pos)
        %    indices(i) = basic_var(q_aux(pos(i)));
        %end
        leaves = min(indices);
    else
        leaves = basic_var(min_pos);
    end
    
    % Reolace the leaving variable for the entering variable in the basis
    basic_var(basic_var == leaves) = enters;
    
end

% the objective value of this optimal basic feasible solution (if the problem
% is bounded)
oval = z_0; 

% vector of size n which is the optimal basic feasible solution 
% corresponding to the obasis.
obfs = zeros(1,n); 
obfs(basic_var) = p;

% vector of size m of indices of column vectors which gives an optimal feasible 
% basis for the problem if the problem is bounded.
obasis = basic_var; 

end
