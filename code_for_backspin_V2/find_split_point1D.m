function [bp1,split_score] = find_split_point1D(sorted_data)

if length(sorted_data(:,1))==1% if dimension is 1 calculate with euclidian distance
    Rcells = squareform(pdist(sorted_data'));
    Rcells = 1- Rcells/std(Rcells(:));
else
    % calc correlation matrix for cells
    Rcells = corr_mat(sorted_data);
end

% look for the optimal breaking point
[N,~] = size(Rcells);
sc = zeros(N,1);
for i=1:N-2
    i
    %     tmp1 = Rcells(1:i,1:i);
    %     tmp2 = Rcells(i+1:N,i+1:N);
    %     sc(i) = (sum(tmp1(:))+sum(tmp2(:)))/ (i^2+(N-i)^2) ;
    if i==1
        tmp1 = Rcells(1:i,1:i);
        tmp1 = sum(tmp1(:));
        tmp2 = Rcells(i+1:N,i+1:N);
        tmp2 = sum(tmp2(:));
        sc(i) = (tmp1+tmp2)/ (i^2+(N-i)^2) ;
    else
        tmp1 = tmp1 + sum(Rcells(i,1:i)) + sum(Rcells(1:i,i));
        tmp2 = tmp2 - sum(Rcells(i,i:end)) - sum(Rcells(i:end,i));
        sc(i) = (tmp1+tmp2)/ (i^2+(N-i)^2) ;
    end
end

[split_score,bp1] = max(sc./mean(sc));
