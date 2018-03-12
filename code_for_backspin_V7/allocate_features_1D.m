% data in linear, working only on the cols
function [ind_gr1_gr2] = allocate_features_1D(data,bp1)

[M,N] = size(data);


Ngr1 = bp1;
Ngr2 = N-bp1;






% divide to two groups
gr1 = 1:bp1;
gr2 = bp1+1:N;
% and assign the features into the two groups
[~,flip_flag1] = min([calc_loccenter(data(:,gr1),1),calc_loccenter(fliplr(data(:,gr1)),1)],[],2);
[~,flip_flag2] = max([calc_loccenter(data(:,gr2),1),calc_loccenter(fliplr(data(:,gr2)),1)],[],2);
data_tmp = data;
data_tmp(flip_flag1==2,gr1) = fliplr(data(flip_flag1==2,gr1));
data_tmp(flip_flag2==2,gr2) = fliplr(data(flip_flag2==2,gr2));

% lc_gr1 = calc_loccenter(data_tmp(:,gr1),1);
% lc_gr2 = calc_loccenter(data_tmp(:,gr2),1);
if Ngr1<=Ngr2    
    loc_center = calc_loccenter([data_tmp(:,gr1),zeros(M,Ngr2-Ngr1),data_tmp(:,gr2)],1);
elseif Ngr1>Ngr2 
    loc_center = calc_loccenter([data_tmp(:,gr1),zeros(M,Ngr1-Ngr2),data_tmp(:,gr2)],1);
end

% loc_center = calc_loccenter(data_tmp,1);
% ind_gr1_gr2 = zeros(size(loc_center));
% ind_gr1_gr2(loc_center<=bp1 ) = 1;
% ind_gr1_gr2(loc_center>bp1) = 2;

ind_gr1_gr2 = zeros(size(loc_center));
ind_gr1_gr2(loc_center<=max([Ngr1,Ngr2]) ) = 1;
ind_gr1_gr2(loc_center>max([Ngr1,Ngr2])) = 2;
if sum(ind_gr1_gr2==1)==0
    d1 = mean(data(:,gr1),2)-mean(data(:,gr2),2);
    [~,imax] = max(d1);
    ind_gr1_gr2(imax) = 1;
elseif sum(ind_gr1_gr2==2)==0
    d1 = mean(data(:,gr2),2)-mean(data(:,gr1),2);
    [~,imax] = max(d1);
    ind_gr1_gr2(imax) = 2;
end
    
    
% % % % % % % % % % % % % % % % % % % % % % % % %