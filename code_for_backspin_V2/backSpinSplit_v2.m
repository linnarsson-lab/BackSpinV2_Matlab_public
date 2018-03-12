% data should be in molecules counts

function [dataout_sorted,genes_order,cells_order,genes_gr_level,cells_gr_level,genes_bor_level,cells_bor_level] = ...
    backSpinSplit_v2(data,numLevels,Nfeture,Nfeture1,N_to_backspin,N_to_cells,mean_tresh,fdr_th,min_gr_cells,min_gr_genes)

[N,M] = size(data);%dimensions of input matrix data
% initialize variables
genes_bor_level = cell(numLevels,1);
cells_bor_level = cell(numLevels,1);

genes_gr_level = zeros(N,numLevels+1);
cells_gr_level = ones(M,numLevels+1);

Nfeture = [Nfeture1,repmat(Nfeture,1,10)];%round(500*1.5.^-[0:10]);%
% N_to_backspin = 500;
% N_to_cells = 500;
% mean_tresh = 0.1;
% fdr_th = 0.1;
% min_gr_cells = 10;
% min_gr_genes = 10;

for i=1:numLevels;% loop over the number of splits requested
    k=0;
    allow_genes_all = find(sum(genes_gr_level(:,i),2)==0);%these are the genes that were not allocated
    log2_m_mat = zeros(length(allow_genes_all),length(unique(cells_gr_level(:,i))));
%     log2_cv_mat = zeros(length(allow_genes_all),length(unique(cells_gr_level(:,i))));
    if i>1
        for j=1:length(unique(cells_gr_level(:,i)));% this loop just calculate the mean expression over the different groups to distribute the genes correctly
            cells_curr = find(cells_gr_level(:,i)==j);
            m_v = mean(data(allow_genes_all,cells_curr),2);
%             cv_v = std(data(allow_genes_all,cells_curr),[],2)./m_v;
            log2_m_mat(:,j) = log2(m_v);
%             log2_cv_mat(:,j) = log2(cv_v);
            log2_m_mat(m_v==0,j) = min(log2_m_mat(~isinf(log2_m_mat(:,j)),j));
%             log2_cv_mat(m_v==0,j) = min(log2_cv_mat(:,j));
        end
        [~,max_m_in] = max(log2_m_mat,[],2);% maximum over the groups
        putative_gene = zeros(size(genes_gr_level(:,i)));
        putative_gene(allow_genes_all) = max_m_in;
    else% in the first round all genes are putative
        putative_gene = ones(N,1);
    end
    for j=1:length(unique(cells_gr_level(:,i))); % loop over the groups of the current splits level
        cells_curr = find(cells_gr_level(:,i)==j);% first take the cell indices
        allow_genes = find( (putative_gene==j | genes_gr_level(:,i)==j) & mean(data(:,cells_curr),2)>mean_tresh);% find the set of allowed genes
        if (length(cells_curr)<min_gr_cells | length(allow_genes)<min_gr_genes) & i>1 ;% skip the split if the groups are too small
            genes_gr_level(genes_gr_level(:,i)==j,i+1) = k+1;
            cells_gr_level(cells_curr,i+1) = k+1;
            k = k+1;
        else
            if length(cells_curr)<=N_to_backspin % if it is normal backspin the genes are just the genes from previous cycle
                genes_curr = find(genes_gr_level(:,i)==j);
            else % otherwise feature selection is used
                m_v = mean(data(allow_genes,cells_curr),2);
                cv_v = std(data(allow_genes,cells_curr),[],2)./m_v;
                log2_m = log2(m_v);
                log2_cv = log2(cv_v);
                x0 = [-0.5,1];
                [param_fit,fval] =  run_min_cvfit(log2_m,log2_cv,x0); 
                param_fit = round(param_fit*100)/100;
                log2_cv_fit = log2( (2.^log2_m).^param_fit(1) + param_fit(2));
                tmp = log2_cv - log2_cv_fit;
                [~,xi1] = sort(tmp,'descend');
                corr_filt = xi1(1:min([Nfeture(i),length(xi1)]));% select according to the number of features setermine by the user
                genes_curr = allow_genes(corr_filt);%[find(genes_gr_level(:,i)==j);allow_genes(corr_filt)];%
            end
            data_tmp = data(genes_curr,cells_curr);% data of the current block (cells x genes)
            data_tmp_log = log2_center(data_tmp);
            if length(cells_curr)>N_to_cells % in this case spilts are done on genes
                [~,genesorder] = sort_and_noplot(data_tmp_log,1,2,5,0);%sorting by spin only genes
                data_tmp_log = data_tmp_log(genesorder,:);%reorder
                data_tmp = data_tmp(genesorder,:);
                genes_curr_sorted = genes_curr(genesorder);
                [bp1,split_score] = find_split_point1D(data_tmp_log');%find break point for 1D order of the genes
                if split_score<0
                    bp1 = 1;
                end
                gr1 = 1:bp1;
                gr2 = bp1+1:length(genes_curr_sorted);
                if ~(bp1==1 | bp1==length(data_tmp(:,1)))
                    data_tmp_log_cn = cent_norm(data_tmp_log); % cent norm the data
                    [~,ind_gr1_gr2] = max([mean(data_tmp_log_cn(gr1,:));mean(data_tmp_log_cn(gr2,:))]);% allocate the cells by the maximum mean over the two groups
                    % % % % % % % % % %
                    if length(cells_curr)>N_to_backspin
                        [~,p] = ttest2(data_tmp_log(:,ind_gr1_gr2==1)',data_tmp_log(:,ind_gr1_gr2==2)');% t test between the groups
                        in = fdr_proc(p,fdr_th);% find the indices that pass FDR
                        genes_curr_sorted = genes_curr_sorted(in);
                        [~,genesorder] = sort_and_noplot(data_tmp_log(in,:),1,2,5,0);%sorting by spin only genes
                        data_tmp_log = data_tmp_log(in(genesorder),:);%reorder
                        data_tmp = data_tmp(in(genesorder),:);
                        genes_curr_sorted = genes_curr_sorted(genesorder);
                        [bp1,split_score] = find_split_point1D(data_tmp_log');%find break point for 1D order of the genes
                        gr1 = 1:bp1;
                        gr2 = bp1+1:length(genes_curr_sorted);
                        data_tmp_log_cn = cent_norm(data_tmp_log);
                        [~,ind_gr1_gr2] = max([mean(data_tmp_log_cn(gr1,:),1);mean(data_tmp_log_cn(gr2,:),1)]);% again, allocate the cells by the maximum mean over the two groups
                        % % % % % % % kNN slassifier
                        numneibor = min(max(ceil(0.01*min([sum(ind_gr1_gr2==1),sum(ind_gr1_gr2==2)])),1),10);% sharpen the group allocation with the KNN two groups classifier
                        ind_gr1_gr2 = knn_corr_classi(data_tmp_log,ind_gr1_gr2,numneibor);   
                        % % % % % % %
                    end
                    % % % % % % % % % % regeister the two groups fro the
                    % next split
                    genes_gr_level(genes_curr_sorted(1:bp1),i+1) = k+1;
                    genes_gr_level(genes_curr_sorted(bp1+1:end),i+1) = k+2;
                    cells_gr_level(cells_curr(ind_gr1_gr2==1),i+1) = k+1;
                    cells_gr_level(cells_curr(ind_gr1_gr2==2),i+1) = k+2;
                    k = k+2;
                else % if the split poin was not valid skip without splitting
                    genes_gr_level(genes_curr_sorted,i+1) = k+1;
                    cells_gr_level(cells_curr,i+1) = k+1;
                    k = k+1;
                end
            else % in this case spilts are done on cells
                [~,~,cellorder] = sort_and_noplot(data_tmp_log,1,3,5,0);%sorting by spin only cells
                data_tmp_log = data_tmp_log(:,cellorder);%reorder
                data_tmp = data_tmp(:,cellorder);
                cells_curr_sorted = cells_curr(cellorder);
                [bp1,split_score] = find_split_point1D(data_tmp_log);%find break point for 1D order of the cells
                % % % % %                 
                ind_gr1_gr2_cells = ones(length(data_tmp_log(1,:)),1);
                ind_gr1_gr2_cells(bp1+1:end) = 2;
                % % % % %                 
                if length(cells_curr)>N_to_backspin
                    [~,p] = ttest2(data_tmp_log(:,1:bp1)',data_tmp_log(:,bp1+1:end)');% t test between the groups
                    in = fdr_proc(p,fdr_th);% find the indices that pass FDR
                    if ~isempty(in)
                        genes_curr = genes_curr(in);
                    else
                        [~,in] = sort(p);
                        genes_curr = genes_curr(in(1:2));
                    end
                    data_tmp = data(genes_curr,cells_curr);
                    data_tmp_log = log2_center(data_tmp);
                    [~,~,cellorder] = sort_and_noplot(data_tmp_log,1,3,5,0);%sorting by spin only cells
                    data_tmp_log = data_tmp_log(:,cellorder);
                    data_tmp = data_tmp(:,cellorder);
                    cells_curr_sorted = cells_curr(cellorder);
                    [bp1,split_score] = find_split_point1D(data_tmp_log);%find break point for 1D order of the cells
                    ind_gr1_gr2_cells = ones(length(data_tmp_log(1,:)),1);
                    ind_gr1_gr2_cells(bp1+1:end) = 2;
                    % % % % % % % kNN slassifier
                    numneibor = min(max(ceil(0.01*min([sum(ind_gr1_gr2_cells==1),sum(ind_gr1_gr2_cells==2)])),1),10);
                    ind_gr1_gr2_cells = knn_corr_classi(data_tmp_log,ind_gr1_gr2_cells,numneibor);
                    % % % % % % %
                end                
                if split_score<0
                    bp1 = 1;
                end
                ind_gr1_gr2 = allocate_features_1D(data_tmp,bp1);%allocate the genes to the two groups
                
                if ~(bp1==1 | bp1==length(data_tmp(1,:)) | sum(ind_gr1_gr2_cells==1)==0 | sum(ind_gr1_gr2_cells==2)==0 | sum(ind_gr1_gr2==1)==0 | sum(ind_gr1_gr2==2)==0)%check for a valid split point
                    if length(cells_curr)<=N_to_backspin
                        cells_gr_level(cells_curr_sorted(1:bp1),i+1) = k+1;
                        cells_gr_level(cells_curr_sorted(bp1+1:end),i+1) = k+2;
                    else
                        cells_gr_level(cells_curr_sorted(ind_gr1_gr2_cells==1),i+1) = k+1;
                        cells_gr_level(cells_curr_sorted(ind_gr1_gr2_cells==2),i+1) = k+2;
                    end
                    genes_gr_level(genes_curr(ind_gr1_gr2==1),i+1) = k+1;
                    genes_gr_level(genes_curr(ind_gr1_gr2==2),i+1) = k+2;
                    k = k+2;
                else % otherwise skip
                    genes_gr_level(genes_curr,i+1) = k+1;
                    cells_gr_level(cells_curr_sorted,i+1) = k+1;
                    k = k+1;
                end
            end
        end
    end
    
end

inc_genes = find(genes_gr_level(:,numLevels+1)>0); % genes that are included are those took part in the last split
%  get a rough order and reorder
[~, xicells] = sortrows(cells_gr_level,[(numLevels+1):-1:1]);
[~, xigenes] = sortrows(genes_gr_level(inc_genes,:),[(numLevels+1):-1:1]);
cells_gr_level = cells_gr_level(xicells,:);
genes_gr_level = genes_gr_level(inc_genes(xigenes),:);

genes_order = inc_genes(xigenes);
cells_order = xicells;
dataout_sorted = data(genes_order,cells_order);
% % % % % % % % % % % % % % % % % % % % % % % % % % order the blocks such
% that similar are closer
gr_uni = unique(cells_gr_level(:,numLevels+1));
num_gr = length(gr_uni);
med_exp = zeros(length(genes_order),num_gr);
for i=1:num_gr % calculate the median expression for each group
    med_exp(:,i) = median(dataout_sorted(:,cells_gr_level(:,numLevels+1)==gr_uni(i)),2);
end
[~,~,gr_cellorder] = sort_and_noplot(log2_center(med_exp),1,3,5,0); % get a 1D order for the groups as blocks
k1 = 0;
k2 = 0;
cellorder = zeros(size(cells_order));
geneorder = zeros(size(genes_order));
for i=1:num_gr % order into blocks
    gr_uni(gr_cellorder(i))
    g_ind = find(genes_gr_level(:,numLevels+1)==gr_uni(gr_cellorder(i)));
    c_ind = find(cells_gr_level(:,numLevels+1)==gr_uni(gr_cellorder(i)));
    cellorder(k1+1:k1+length(c_ind)) = c_ind;
    geneorder(k2+1:k2+length(g_ind)) = g_ind;
    k1 = k1+length(c_ind);
    k2 = k2+length(g_ind);
end
geneorder(geneorder==0) =[];
cellorder(cellorder==0) =[];
genes_order = genes_order(geneorder);
cells_order = cells_order(cellorder);
dataout_sorted = dataout_sorted(geneorder,cellorder);
cells_gr_level = cells_gr_level(cellorder,:);
genes_gr_level = genes_gr_level(geneorder,:);
% % % % % % % % % % % % % % % % % % % % % % % % % %resort each block
% finally run SPIN on each block the get the local order of genes and cells
% (within the block)
for i=1:num_gr
    g_ind = find(genes_gr_level(:,numLevels+1)==gr_uni(i));
    c_ind = find(cells_gr_level(:,numLevels+1)==gr_uni(i));
    if length(g_ind)>1 & length(c_ind)>1
        data_tmp_log = log2_center(dataout_sorted(g_ind,c_ind));
        [~,genesorder,cellorder] = sort_and_noplot(data_tmp_log,1,1,5,0);
        dataout_sorted(g_ind,:) = dataout_sorted(g_ind(genesorder),:);
        dataout_sorted(:,c_ind) = dataout_sorted(:,c_ind(cellorder));
        genes_order(g_ind,:) = genes_order(g_ind(genesorder),:);
        cells_order(c_ind,:) = cells_order(c_ind(cellorder),:);
        genes_gr_level(g_ind,:) = genes_gr_level(g_ind(genesorder),:);
        cells_gr_level(c_ind,:) = cells_gr_level(c_ind(cellorder),:);
    end
end



for i=1:numLevels
    genes_bor_level{i} = find(diff(genes_gr_level(:,i+1))~=0)+1;
    cells_bor_level{i} = find(diff(cells_gr_level(:,i+1))~=0)+1;
end






































