


function [x,fval] =  run_min_cvfit(log2_m,log2_cv,x0) 
log2_m1 = log2_m(:);
log2_cv1 = log2_cv(:);
if length(log2_m1)>1000    
    [n,xi] = hist(log2_m1,50);
    xi = [min(log2_m1),xi,max(log2_m1)];
    med_n = prctile(n,50);
%     med_n = prctile(n,20);
    for i=2:(length(n)+1)
        ind = find(log2_m1>xi(i-1) & log2_m1<xi(i));
        if length(ind)>med_n
            rng('default');
            ind = ind(randperm(length(ind)));
            ind = ind(1:round(length(ind)-med_n));
            log2_m1(ind) = [];
            log2_cv1(ind) = [];
        elseif round(med_n/length(ind))>1 & length(ind)>3
            log2_m1 = [log2_m1;repmat(log2_m1(ind),round(med_n/length(ind))-1,1)];
            log2_cv1 = [log2_cv1;repmat(log2_cv1(ind),round(med_n/length(ind))-1,1)];
        end
    end
end
            
            

[x,fval] = fminsearch(@nestedfun,x0,optimset('MaxIter',1e20,'TolFun',1e-9,'tolx',1e-9,'display','off'));
% Nested function that computes the objective function     
    function y = nestedfun(x)
        alpha = x(1);
        k = x(2);
%         y = sum((log2( (2.^log2_m1).^alpha +k) - log2_cv1).^2);
        y = sum(abs((log2( (2.^log2_m1).^alpha +k) - log2_cv1)));
%         y = sum(abs((log2( (2.^log2_m1).^alpha)+k - log2_cv1)));
    end
end








