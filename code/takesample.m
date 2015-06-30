function takesample(k)
    n = 100;
    T = nchoosek(1:n, k);
    Z = nchoosek(n,k);
%     T_file = ['T-',num2str(k)];
%     dlmwrite(T_file, [k*ones(Z,1),T], '\t');
%     dlmwrite([T_file, '.stat'], [n; 1], '\t');

    lens = [1,2,5,10,20,50,100,200,500,1000,2000];
    

    for l=lens
        for idx=0:9
            p = 1.0/Z;
            L = binornd(l*ones(Z,1), p);

            S = [];
            for i=1:Z
                S = [S; repmat(T(i,:), L(i),1)];
            end

            S_file = ['S-', num2str(k), '-', num2str(l), '-', num2str(idx)];
            dlmwrite(S_file, [k*ones(size(S,1),1), S], '\t');
            dlmwrite([S_file, '.stat'], [n; l], '\t');
        end
    end
end