% INPUT:
% ------
%   A       - (n�n double) matrix
%   b       - (n�1 double) vector
%
% -------
% OUTPUT:
% -------
%   x       - (n�1 double) solution of the linear system Ax=b
%
% -----
% NOTE:
% -----
%   --> This function is meant to illustrate Gaussian elimination with
%       partial pivoting.
%
%==========================================================================
function x = Gaussian_Elimination(A,b)
    
    % determines n
    n = size(A,1);
    
    % creates augmented matrix
    A = [A,b];
    
    % keeps track if matrix is singular
    singular = false;
    
    % elimination process
    for i = 1:(n-1)
        
        % determines pivot row
        p = i:1:n;
        p(A(p,i)==0) = [];
        p = min(p);
        
        % if all possible pivots in the column are zero the matrix is singular
        if max(abs(A(:,i))) <= eps
            singular = true;
            break
        end
        
        % if p does not equal i, switches the ith and pth rows
        if p ~= i
            Ai = A(i,:);
            Ap = A(p,:);
            A(i,:) = Ap;
            A(p,:) = Ai;
        end
        
        % elementary row operation
        for j = (i+1):n
            A(j,:) = A(j,:)-(A(j,i)/A(i,i))*A(i,:);
        end
        
    end
    
    % if bottom right element of A is 0,  
    % Then A is singular
    if abs(A(n,n)) <= eps
        singular = true;
    end
    
    
    x = zeros(n,1);
    
    % performs backward substitution to solve Ax = b if matrix nonsingular
    if ~singular
        x(n) = A(n,n+1)/A(n,n);
        for i = (n-1):(-1):1
            S = 0;
            
            %vectorize version
            S = sum(A(i,(i+1):n).* transpose(x((i+1):n)));
            
            %for j = (i+1):n
            %    S = S+A(i,j)*x(j);
            %end
            
            x(i) = (A(i,n+1)-S)/A(i,i);
        end
    end
   
    if singular
        warning('Matrix is singular to working precision.')
    end
    
end