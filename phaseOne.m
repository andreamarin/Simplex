% First Assignment.  
%   Andrea Marin Alarcon 158999
%   Andrea Perez Vega 154467
%    
% The objective of the code is to solve (if it is possible) the following
% linear program:
%          maximise c^T x
%           subject to Ax = b, x >= 0, b >=0

% In Phase One, we solve an auxiliary program in order to get a starter basic feasible solution 
% and it's corresponding basis, if they exist
function[nvac, basis, bfs] = phaseOne(A, b, c)
% INPUT:
% A mxn matrix with m <= n and rank of A is m.
% b column vector with m rows.
% c column vector with n rows.
% OUTPUT:
% nvac = 0 if the feasible set is empty, 0 otherwise
% basis = feasible basis for the original LP problem
% bfs = basic feasible solution for the original LP problem

    k = length(b); 
    [~,n] = size(A);

    % the indices corresponding to the correction variables in the auxiliary program
    correction_variables = n+1:n+k; 
    
    % First, we make sure that b >= 0 
    indices = find(b<0); % variable indices where b  < 0
    if isempty(indices) == false
        A(indices,:) = -1*A(indices,:);
        b(indices) = -1*b(indices);
    end
    
    %The new matrix A, adding the correction variables;
    A_aux(:,1:n) = A; 
    A_aux(:, correction_variables) = eye(k);

    % The function to be maximised in the auxiliary program
    % where the original variables' coefficients are zero and the function is given by:
    % -x_(n+1) - x_(n+2) - ... - x_(n+k)
    c_aux = zeros(n+k,1); 
    c_aux(correction_variables) = -1 * ones(k,1);

    % The starting basic feasible solution for the auxiliary program
    sbfs = zeros(n+k, 1);
    sbfs(correction_variables) = b;

    %We apply phaseTwo to see if the auxiliary program be solved.
    [bound, obasis, obfs, oval] = phaseTwo(A_aux, b, c_aux, correction_variables, sbfs);

    if(bound == 0)
        nvac = 0; % The feasible set is empty.
    else
        %If the optimum value is zero in the auxiliary program, the feasible set is non-empty.
        if (oval == 0) 
            nvac = 1; %The feasible set is non-empty.
        else
            nvac = 0;
        end
    end
    
    if nvac == 1 && max(obasis) > n
    % If there is a basic feasible solution, we check if we have a degenerate basis 
    % and if so, we build a feasible basis with the original variables
        non_basic = setdiff(1:n,obasis);
        basis = obasis;
        
        i = 1;
        max_index = max(basis);
        while max_index > n 
            
            index = basis == max_index;
            % we replace the correction variable in the basis with a non basic variable
            basis(index) = non_basic(i);
            
            aux_basis = setdiff(basis, correction_variables);
            
            % check if the columns of matrix A corresponding to the basis
            % are linearly independent
            while rank(A(aux_basis)) ~= length(aux_basis)
                i = i+1;
                basis(index) = non_basic(i);
                aux_basis = setdiff(basis, correction_variables);
            end
            
            max_index = max(basis);
        end
        
    else
        basis = obasis;
    end
    
    bfs = obfs(1:n); %Vector of size n of the basic feasible solution corresponding to the original basis.
end

