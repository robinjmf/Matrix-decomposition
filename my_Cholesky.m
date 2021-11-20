function L = my_Cholesky(A)

%Cholesky decomposition is used for positive definite matrices and it can
%factorize The Matrix A to L * L' , 
%this function input matrix as A argument and output the L matrix
%which was defined in previous line 

    if ~is_pos_def(A) == 1
        error('Cholesky Decomposition is used only for positive definite matrices')
    end
    n = rank(A);
    L = zeros(n);
    
    L(1,1) = sqrt(A(1,1));
    L(2:n,1) = A(2:n,1)/L(1,1);
        
    for i = 2:n
            
        L(i,i)= sqrt(A(i,i)- sum(power(L(i,1:i-1),2)));
        %vectorize version
        L(i+1:n,i) = (A(i,i+1:n).' - (L(i+1:n,1:i-1) * L(i,1:i-1).')) / L(i,i);    
       
        %for j = i+1:n
            %L(j,i) = (A(i,j)-(dot(L(j,1:i-1), L(i,1:i-1))))/L(i,i);
        %end
            
        
    end
end

    
    
    