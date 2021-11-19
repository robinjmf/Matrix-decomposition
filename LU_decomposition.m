function LU = LU_decomposition( A )

%LU decomposition of A: 
%L(i,j) = LU(i,j) --> i>j
%L(i,j) = 0       --> j>i
%L(i,j) = 1       --> i=j
%U(i,j) = LU(i,j) --> j>=i
%U(i,j) = 0       --> i>j

    n = rank(A);
    LU = zeros(n);
    LU(1,1:n) = A(1,1:n);
    LU(2:n,1) = A(1,2:n)./A(1,1);
    
    for k = 2:n
        for j = k:n
            
            LU(k,j) = A(k,j) - dot(LU(k,1:k-1) , LU(1:k-1,j));
        end
        
        for i = k+1:n
            
            LU(i,k) = (A(i,k) - dot(LU(i, 1:k-1) , LU(1:k-1, k)))/LU(k,k);
        end
    end   
end

