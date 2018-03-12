% % lin_log_flag=1 for linear, and lin_log_flag=2 for log2
function [loc_center] = calc_loccenter(x,lin_log_flag)

[M,N] = size(x);
if N==1 & M>1
    x = x';
end
[M,N] = size(x);
loc_center = zeros(M,1);
min_x = min(x,[],2);
x = x - repmat(min_x,1,N);
for i=1:M
    ind = find(x(i,:)>0);
    if ~isempty(ind)
        if lin_log_flag==1
            w = x(i,ind)/sum(x(i,ind));
        else
            w = (2.^x(i,ind))/sum(2.^x(i,ind));
        end
        loc_center(i) = sum(w.*ind);%mean(ind,2);        
    else
        loc_center(i) = 0;
    end
end
