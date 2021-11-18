function f = is_pos_def(input)
    eig_A = eig(input);
    flag = 0;
    for i = 1:rank(input)
        if eig_A(i) <= 0
            flag = 1;
        end
    end
    f = flag;
    if flag == 1
        disp('the matrix is not positive definite')
    else
        disp('the matrix is positive definite')
    end


