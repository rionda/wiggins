function opt(input_sample, number_of_nodes)
tic
    round = 0;
    theta = 0.5;
    d = importdata(input_sample);
    d = d + 1;
    r = size(d,1);
    disp(r)
    S = cell(r,1);
    T = sparse(r, number_of_nodes);
%     N2S=sparse(number_of_nodes, r);
    for i=1:r
        tmp = d(i,:);
        S{i} = tmp(~isnan(tmp));
        T(i, S{i}) = 1;
    end

    function [f, gf] = lbcost(p)
        round = round + 1
%         gf = zeros(number_of_nodes,1);
        tmp = 1.0 ./ (1-theta*(1-T*p));
        f = sum(tmp) / theta;
        gf= T'*(-tmp.^2);
        
%         for j=1:r
%             par = 1.0/(1-theta*())
% %             par = 1.0/(1-theta*(1-sum(p(S{j}))));
%             f = f + par;
%             for k=S{j}
%                gf(k) = gf(k) - par^2;
%             end
%            
%         end
        f = f/theta;
    end
    p = ones(number_of_nodes,1)/number_of_nodes;
    
    
    
options = optimset('Algorithm', 'interior-point');
options = optimset(options, 'GradObj', 'on');
options = optimset(options, 'Display', 'off');
options = optimset(options,'UseParallel','always');


%{
options = optimoptions('fmincon','GradObj','on');
options = optimoptions('fmincon','Hessian','user-supplied');
%}

N = number_of_nodes;
p0 = ones(N,1)/N;

[p, fval] = fmincon(@lbcost, p0, [], [], ones(1,N), 1, zeros(1,N) , ones(1,N), [], options);

fval
toc
% dlmwrite(lowerBound, fval);
% dlmwrite(outputFile, p/sum(p));


% 
% %%%%%%%%%%%%%%%%%%%%%%% shoul I make S sparse!?!
% 
% % Defining my cost function
% function [f, gf] = lbcost(p)
%     
%     denoms_1 = 1./p;
%     denoms_2 = 1 ./ (S * p);
% %     disp(class(rho))
% %     disp(size(myPi))
% %     disp(size(denoms_1))
% %     disp(size(denoms_2))
%     myvec_1 = (1-rho) .* (myPi .* denoms_1);
%     myvec_2 = (rho)   .* (myPi .* denoms_2);
%     f = sum(myvec_1 + myvec_2);
% 
%     gf_1 = - (myvec_1 .* denoms_1);
%     gf_2 = - S * (myvec_2 .* denoms_2);
%     gf = gf_1 + gf_2;
% end
% 
% 
%  
% 
% 
% 
% % then running the optimization
% 
% options = optimset('Algorithm', 'interior-point');
% options = optimset(options, 'GradObj', 'on');
% options = optimset(options, 'Display', 'off');
% options = optimset(options,'UseParallel','always');
% 
% 
% %{
% options = optimoptions('fmincon','GradObj','on');
% options = optimoptions('fmincon','Hessian','user-supplied');
% %}
% 
% p0 = ones(N,1)/N;
% 
% [p, fval] = fmincon(@lbcost, p0, [], [], ones(1,N), 1, zeros(1,N) , ones(1,N), [], options);
% 
% dlmwrite(lowerBound, fval);
% dlmwrite(outputFile, p/sum(p));

end