function LLT = my_Cholesky(A)

    if is_pos_def(A) == 1
        disp('Cholesky Decomposition is used only for positive definite matrixes')
    else
        n = rank(A);
        L = zeros(rank(A));
        L(1,1) = sqrt(A(1,1));
        
        for i=2:n
            L(i,1)=A(i,1)/L(1,1);
        end
        
        for i = 2:n
            
            L(i,i)= sqrt(A(i,i)- sum(power(L(i,1:i-1),2)));
            for j = i+1:n
                L(j,i) = (A(i,j)-(dot(L(j,1:i-1), L(i,1:i-1))))/L(i,i);
            end
            
        end
    end
    LLT = L;
end

    
    
    