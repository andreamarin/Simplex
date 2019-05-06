% First Assignment.  
%   Andrea Marin Alarcon 158999
%   Andrea Perez Vega 154467
%    
% The objective of the code is to solve (if it is possible) the following
% linear program:
%          maximise c^T x
%           subject to Ax = b, x >= 0, b >=0

% In this part, we implement the entire simplex algorithm using both of the previous phases, this function
% tells us if the problem has a basic feasible solution, if it is unbounded,
% and finds the optimal solution if it exists.
function[status, obasis, obfs, oval] = bothPhases(A, b, c)
% INPUT:
% A mxn matrix with m <= n and rank of A is m.
% b column vector with m rows.
% c column vector with n rows.
% OUTPUT:
% status = -1 if the feasible set is empty
% status = 0 if the feasible set is non-empty but the problem is unbounded
% status = 1 if the problem is bounded (there is an optimal solution)
% obasis = a vector of size m of indices of an optimal feasible basis for the problem 
% obfs = a vector of size n which is the optimal basic feasible solution 
% corresponding to this optimal basis if the feasible set is non-empty and the problem is bounded
% oval = the objective value of this optimal basic feasible solution 

% First we find a basic feasible solution for the problem
[nvac, sbasis, sbfs] = phaseOne(A, b, c);

if (nvac == 0)
    status = -1; %Feasible set is empty, then it can't be solved.
    obasis = sbasis;
    obfs = sbfs;
    oval = -1;
else
    [bound, obasis, obfs, oval] = phaseTwo(A, b, c, sbasis, sbfs);
    % obasis: Vector of size m of indices of an optimal feasible basis for the problem if the feasible set is 
    % non-empty and the problem is bounded.
    % obfs: Vector of size n which is the optimal basic feasible solution corresponding to the obasis. 
    % oval: The objective value of this obfs.
    
    if(bound == 0)
        status = 0; %The feasible set is non-empty but the problem is ubounded (not optimal solution).
    else
        status = 1; %The problem set is bounded, therefore, has an optimal solution.
    end    
end

end
   