%Andrea Marin Alarc�n 158999
%Andrea P�rez Vega 154467

function[bound, obasis, obfs, oval] = firstassignment(A, b, c, sbasis, sbfs)
%maximise c^T x
% subject to Ax = b, x >= 0, b >=0
%
% Input:
% A mxn matrix with m <= n and rank of A is m
% b column vector with m rows
% c column vector with n rows
% sbasis a vector of size m of indices of column vectors for a feasible basis for this
%problem from which to start the simplex method
% sbfs a vector of size n which is the basic feasible solution corresponding to this
%basis
%
% Output:
% bound = 0 if the problem is unbounded (there is no optimal solution)
% bound = 1 if the problem is bounded (there is an optimal solution)

% obasis = a vector of size m of indices of column vectors which gives an optimal fea-
%sible basis for the problem if the problem is bounded

% obfs = a vector of size n which is the optimal basic feasible solution correspond-
%ing to this optimal basis if the problem is bounded

% oval = the objective value of this optimal basic feasible solution (if the problem
%is bounded)

[m,n] = size(A);

A_b = zeros(m,m);
A_n = zeros(m, n-m);
c_b = zeros(m,1);
c_n = zeros(n-m,1);

finished = false;

basic_var = sbasis;

while finished == false
    non_basic = setdiff(1:n, basic_var);

    for i = 1:1:m
        index = basic_var(i);
        A_b(:,i) = A(:,index);
        c_b(i) = c(index);
    end

    for i = 1:1:(n-m)
        index = non_basic(i);
        A_n(:,i) = A(:,index);
        c_n(i) = c(index);
    end

    Ab_inverse = inv(A_b);

    Q = -Ab_inverse * A_n;
    p = Ab_inverse * b;
    z_0 = transpose(c_b)* Ab_inverse * b;
    r = c_n - transpose(transpose(c_b) * Ab_inverse * A_n)

    %We will use Bland's rule to avoid cycles
    r_aux = find(r > 0);
    
    if isempty(r_aux)
        finished = true;
        bound = 1;
        continue
    end
    
    indices = zeros(1, length(r_aux));
    for i = 1:1:length(r_aux)
        index = r_aux(i);
        indices(i) = non_basic(index);
    end
    enters = min(indices)

    q_aux = find(Q(:,enters) < 0);
    if isempty(q_aux)
        finished = true;
        bound = 0;
        continue
    end
    
    constraints = zeros(1, length(q_aux))
    for i = 1:1:length(q_aux)
        index = q_aux(i)
        constraints(i) = -p(index)/Q(index,enters)
    end

    min_constraint = min(constraints);
    pos = find(constraints == min_constraint);
    
    if length(pos) > 1
        indices = zeros(length(pos),1)
        for i = 1:1:length(pos)
            indices(i) = basic_var(pos(i))
        end
        leaves = min(indices);
    else
        leaves = basic_var(pos(1));
    end
    
    index = find(basic_var == leaves);
    basic_var(index) = enters
    
end

oval = z_0;
obfs = zeros(1,n);
obasis = basic_var;
for i = 1:1:m
    index = basic_var(i);
    obfs(index) = p(i);
end

end
    


    
