function [bp1,split_score] = find_split_point1D_v4(sorted_data,stop_th,min_gr_cells,flag_val_stip)

if length(sorted_data(:,1))==1% if dimension is 1 calculate with euclidian distance
    Rcells = squareform(pdist(sorted_data'));
    Rcells = 1- Rcells/std(Rcells(:));
else
    % calc correlation matrix for cells
    Rcells = corr_mat(sorted_data);
end
Rcells1 = Rcells;
Rcells(Rcells<0) = 0;
[N,M] = size(sorted_data);

N = length(Rcells);
if N<30
    stop_th1 = stop_th(1);
elseif N>200
    stop_th1 = stop_th(2)
else
    stop_th1 = interp1([30,200],stop_th,N);
end
wid = 1;
w_pow = 2;
[I,J] = meshgrid(1:N,1:N);
weights_mat = abs(exp(-((I-J).^w_pow)/wid/N)) - 1e-8;%elimnate values that are too small and make copmutation slow
weights_mat(weights_mat<0) = 0;
wid = wid/(mean(weights_mat(:))/mean(Rcells(:)))^2;
wid(wid<0.01) = 0.01;
weights_mat = abs(exp(-((I-J).^w_pow)/wid/N)) - 1e-8;
weights_mat(weights_mat<0) = 0;
weights_mat = weights_mat/mean(weights_mat(:))*mean(Rcells(:));


if isempty(Rcells)
    bp1 = [];
    split_score = 0;
else
    % look for the optimal breaking point
    [N,~] = size(Rcells);
    sc = zeros(N-2,1);
    for i=1:N-2
        i
        if i==1
            tmp1 = Rcells(1:i,1:i);
            tmp1 = sum(tmp1(:))-i;
            tmp2 = Rcells(i+1:N,i+1:N);
            tmp2 = sum(tmp2(:))-(N-i);
            sc(i) = (tmp1+tmp2)/ (-N + i^2+(N-i)^2) ;
        else
            tmp1 = tmp1 + sum(Rcells(i,1:i-1)) + sum(Rcells(1:i-1,i));
            tmp2 = tmp2 - sum(Rcells(i,i+1:end)) - sum(Rcells(i+1:end,i));
            sc(i) = (tmp1+tmp2)/ (-N + i^2+(N-i)^2)  ;
        end
    end
    
    sc_w = zeros(N-2,1);
    for i=1:N-2
        i;
        if i==1
            tmp1 = weights_mat(1:i,1:i);
            tmp1 = sum(tmp1(:))-i;
            tmp2 = weights_mat(i+1:N,i+1:N);
            tmp2 = sum(tmp2(:))-(N-i);
            sc_w(i) = (tmp1+tmp2)/ (-N + i^2+(N-i)^2) ;
        else
            tmp1 = tmp1 + sum(weights_mat(i,1:i-1)) + sum(weights_mat(1:i-1,i));
            tmp2 = tmp2 - sum(weights_mat(i,i+1:end)) - sum(weights_mat(i+1:end,i));
            sc_w(i) = (tmp1+tmp2)/ (-N + i^2+(N-i)^2)  ;
        end
    end
    sc_w_smooth = smoothn(sc_w,max(min(ceil(length(sc_w)*0.05),10),3));
    sc_smooth = smoothn(sc,max(min(ceil(length(sc)*0.1),10),3));
    line1 = sc./sc_w;
    line1_smooth = sc_smooth./sc_w_smooth;
    d1_line1 = d1calc([1:length(line1)]',line1);
    d1_line1_smooth = d1calc([1:length(line1)]',line1_smooth);
    d2_line1 = d2calc([1:length(line1)]',line1);
    d2_line1_smooth = d2calc([1:length(line1)]',line1_smooth);
    %     maxpoints = find(diff(sign(diff(line1)))<0) + 1;
    %     maxpoints_smooth = (find(diff(sign(diff(line1_smooth)))<0)) +1;
    maxpoints = find(diff(sign(d1_line1))<0);
    maxpoints_smooth = (find(diff(sign(d1_line1_smooth))<0));
    if N>100
        maxpoints_smooth(maxpoints_smooth<N*0.01 | maxpoints_smooth>N*0.99) = [];
    end
    maxpoints_smooth(maxpoints_smooth<min_gr_cells | maxpoints_smooth>(N-min_gr_cells)) = [];
    stipness_max_points = d2_line1_smooth(maxpoints_smooth)
    %     if length(maxpoints_smooth)>1
    %         [zz,stip_xi] = sort(abs(stipness_max_points),'descend');
    %         imax = maxpoints_smooth(stip_xi(1));
    %     else
    %         [~,imax] = max(line1_smooth(maxpoints_smooth));
    %         imax = maxpoints_smooth(imax);
    %     end
    if flag_val_stip==1
        [~,imax] = max(line1_smooth(maxpoints_smooth));
        imax = maxpoints_smooth(imax);
        split_score = line1_smooth(imax)
        bp1 = imax+1
    elseif flag_val_stip==2
        [~,imax] = max(stipness_max_points);
        imax = maxpoints_smooth(imax);
        split_score = line1_smooth(imax)
        bp1 = imax+1
    end
%         stip_max = mean(abs([line1_smooth(imax)-line1_smooth(imax-10),line1_smooth(imax)-line1_smooth(imax+10)]))
    
%     h111 = figure(111); set(gcf,'position',[200,400,600,500]); imagesc(Rcells1); hold on; plot(maxpoints_smooth,maxpoints_smooth,'sk'); plot(bp1,bp1,'pm','markersize',12);
%     h222 = figure(222); set(gcf,'position',[800,400,600,500]); plot(line1); hold on; plot(line1_smooth,'k'); plot(maxpoints_smooth,line1_smooth(maxpoints_smooth),'or');
%     [N,M] = size(sorted_data)
    if split_score<stop_th1 %split_score*sqrt(N)*num_mxi<stop_th ;%split_score<stop_th/num_mxi;%
        bp1 = [];
    end
    if isempty(split_score);
        fprintf('stop');
    end
%     close(h111)
%     close(h222)
end