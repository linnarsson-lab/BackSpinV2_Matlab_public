
tic
clear all
close all



load /data2/C1_stuff/hypothalamus_stress_analysis/afterLoading_analysis_hypoth_fromdatabase
% % % % % % % % % % % % % 
geneuni = genes_hypoth;

total_mol = sum(moldata_hypoth);

uni_platevec = unique(chip_hypoth);
samples_all = cell(length(chip_hypoth),1);
for i=1:length(chip_hypoth)
    samples_all{i} = [chip_hypoth{i},'_',well_hypoth{i}];
end

% % % % % exclude low wells
goodwells = total_mol>1.5e3 & total_mol<70e3 & ~(strcmpi(chip_hypoth,'1772078075') | strcmpi(chip_hypoth,'1772075316'));
moldata_1 = moldata_hypoth(:,goodwells);
chip_1 = chip_hypoth(:,goodwells);
well_1 = well_hypoth(:,goodwells);
age_1 = age_hypoth(:,goodwells);
sex_1 = sex_hypoth(:,goodwells);%% female=1, male =-1, mix=0;
diam_1 = diam_hypoth(:,goodwells);
stress6h_1 = stress6hr_hypoth(:,goodwells);
tissue_all_1 = tissue_hypoth(goodwells);
total_mol_1 = total_mol(:,goodwells);
image_all_1 = image_all_hypoth(goodwells);
moldata_spikes = moldata_hypoth_ercc(:,goodwells);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


tissue_uni = unique(tissue_all_1);

indgood = mean(moldata_1,2)>(50/length(chip_1)) & sum(moldata_1>0,2)<length(chip_1)*0.7 & sum(moldata_1>0,2)>10;
moldata_1 = moldata_1(indgood,:);
geneuni_1 = geneuni(indgood,:);


cellid_1 = cell(length(chip_1),1);
for i=1:length(chip_1)
    cellid_1{i} = [chip_1{i},'_',well_1{i}];
end



load /data2/C1_stuff/hypothalamus_stress_analysis/data_after_backspinv2_lev8_hypoth_04-Jun-2015.mat
% contain variables: 'moldata_1_backspin7', 'cellid_1', 'genes_1','genes_bor_level','cells_bor_level','genes_gr_level',
% 'cells_gr_level','splitlev','Nfeture','N_to_backspin','N_to_cells','mean_tresh','fdr_th','min_gr_cells','min_gr_genes');


cellid_orig = samples_all(goodwells);
[~,loc] = ismember(cellid_1,cellid_orig);

dataout_sorted_all = moldata_1(:,loc);

neurons_ind = [1002:2080];
cell_neurons = cellid_1(neurons_ind);
nonneurons_ind = setdiff([1:length(cellid_1)],neurons_ind);

[~,p] = ttest2(log2(dataout_sorted_all(:,neurons_ind)+1)',log2(dataout_sorted_all(:,nonneurons_ind)+1)','tail','right');
neurons_genes_ind = fdr_proc(p,0.05);
neurons_genes = geneuni_1(neurons_genes_ind);
moldata_neurons = dataout_sorted_all(neurons_genes_ind,neurons_ind);
moldata_neurons_allgenes = dataout_sorted_all(:,neurons_ind);

tic
splitlev = 7;
Nfeture1 = 500;
Nfeture = 200;
N_to_backspin = 10;
N_to_cells = 500;
mean_tresh = 0.1;
fdr_th = 0.3;
min_gr_cells = 5;
min_gr_genes = 3;

[dataout_sorted_n,genes_order_n,cells_order_n,genes_gr_level_n,cells_gr_level_n,genes_bor_level_n,cells_bor_level_n] = ...
    backSpinSplit_v2_fix(moldata_neurons,splitlev,Nfeture,Nfeture1,N_to_backspin,N_to_cells,mean_tresh,fdr_th,min_gr_cells,min_gr_genes);
toc


cell_neurons_sorted = cell_neurons(cells_order_n);
neurons_genes_sorted = neurons_genes(genes_order_n);
% 
% save(['data_after_backspinv2_neurons_lev',num2str(splitlev),'_hypoth_',date],'dataout_sorted_n', 'cell_neurons_sorted', 'neurons_genes_sorted'...
%     ,'genes_bor_level_n','cells_bor_level_n','genes_gr_level_n','cells_gr_level_n','splitlev','Nfeture','N_to_backspin','N_to_cells','mean_tresh','fdr_th','min_gr_cells','min_gr_genes');
load /data2/C1_stuff/hypothalamus_stress_analysis/data_after_backspinv2_neurons_lev7_hypoth_06-Jun-2015.mat

[~,loc] = ismember(cell_neurons_sorted,cell_neurons);
dataout_sorted_all_n = moldata_neurons(:,loc);
sorted_moldata_neurons_allgenes = moldata_neurons_allgenes(:,loc);

sorted_data_cn = cent_norm(log2_center(dataout_sorted_n));

T_cells_tmp = cells_gr_level_n(:,splitlev+1)';
T_genes_tmp = genes_gr_level_n(:,splitlev+1)';

cells_bor_tmpgr = cells_bor_level_n;
genes_bor_tmpgr = genes_bor_level_n;

lev = splitlev;
cells_bor_tmpgr{lev} = [cells_bor_tmpgr{lev};length(T_cells_tmp)+1];
genes_bor_tmpgr{lev} = [genes_bor_tmpgr{lev};length(T_genes_tmp)+1];

figure;
set(gcf,'position',[100,100,600,1000],'color','w')
axes('position',[0.1,0.05,0.85,0.85])
imagesc(sorted_data_cn,[prctile(sorted_data_cn(:),1),prctile(sorted_data_cn(:),99)]);
hold on;
linewid =0.5;
bor_color = 'grey11';%'green1';%
for jj=1:length(cells_bor_tmpgr{lev})
    if jj>1
        plot(cells_bor_tmpgr{lev}(jj)*[1,1]-0.5,[genes_bor_tmpgr{lev}(jj-1),genes_bor_tmpgr{lev}(jj)],'-','linewidth',linewid,'color',get_RGB(bor_color))
        plot([cells_bor_tmpgr{lev}(jj-1),cells_bor_tmpgr{lev}(jj)]-0.5,genes_bor_tmpgr{lev}(jj)*[1,1] ,'-','linewidth',linewid,'color',get_RGB(bor_color))
        plot(cells_bor_tmpgr{lev}(jj-1)*[1,1]-0.5,[genes_bor_tmpgr{lev}(jj-1),genes_bor_tmpgr{lev}(jj)]-0.5,'-','linewidth',linewid,'color',get_RGB(bor_color))
        plot([cells_bor_tmpgr{lev}(jj-1)-0.5,cells_bor_tmpgr{lev}(jj)],genes_bor_tmpgr{lev}(jj-1)*[1,1]-0.5 ,'-','linewidth',linewid,'color',get_RGB(bor_color))
    else
        plot(cells_bor_tmpgr{lev}(jj)*[1,1]-0.5,[1,genes_bor_tmpgr{lev}(jj)]-0.5,'-','linewidth',linewid,'color',get_RGB(bor_color))
        plot([1,cells_bor_tmpgr{lev}(jj)]-0.5,genes_bor_tmpgr{lev}(jj)*[1,1] -0.5,'-','linewidth',linewid,'color',get_RGB(bor_color))
        plot(1*[1,1]-0.5,[1,genes_bor_tmpgr{lev}(jj)]-0.5,'-','linewidth',linewid,'color',get_RGB(bor_color))
        plot([1,cells_bor_tmpgr{lev}(jj)]-0.5,1*[1,1] -0.5,'-','linewidth',linewid,'color',get_RGB(bor_color))
    end
end
set(gca,'xtick',[0:200:3000],'ytick',[0:100:length(T_genes_tmp)])
colormap('summer');freezeColors(gca);

cellid_orig = samples_all(goodwells);
[~,loc] = ismember(cell_neurons_sorted,cellid_orig);

sorted_age_1 = age_1(loc);
sorted_sex_1 = sex_1(loc);
sorted_diam_1 = diam_1(loc);
sorted_total_mol_1 = total_mol_1(loc);
sorted_stress_1 = stress6h_1(loc);
sorted_chip_1 = chip_1(loc);
sorted_chip_1_uni = unique(sorted_chip_1);
sorted_chip_num = zeros(size(chip_1));
for i=1:length(sorted_chip_1_uni)
    sorted_chip_num(strcmpi(sorted_chip_1,sorted_chip_1_uni{i})) = i;
end

sorted_age_1(isnan(sorted_age_1)) = mean(sorted_age_1(~isnan(sorted_age_1)));

axes('position',[0.1,0.9,0.85,0.01])
imagesc(cent_norm(T_cells_tmp))
set(gca,'xtick',[],'ytick',1,'yticklabel','clusters')
colormap('jet');freezeColors(gca);
axes('position',[0.1,0.91,0.85,0.01])
imagesc(cent_norm(sorted_age_1))
set(gca,'xtick',[],'ytick',1,'yticklabel','age')
colormap('jet');freezeColors(gca);
axes('position',[0.1,0.92,0.85,0.01])
imagesc(cent_norm(sorted_chip_num))
set(gca,'xtick',[],'ytick',1,'yticklabel','chip')
colormap('jet');freezeColors(gca);
axes('position',[0.1,0.93,0.85,0.01])
imagesc(-cent_norm(sorted_sex_1))
set(gca,'xtick',[],'ytick',1,'yticklabel','sex')
colormap('gray');freezeColors(gca);

axes('position',[0.1,0.94,0.85,0.01])
imagesc(cent_norm(sorted_diam_1))
set(gca,'xtick',[],'ytick',1,'yticklabel','diam')

axes('position',[0.1,0.95,0.85,0.01])
imagesc(-cent_norm(sorted_stress_1))
set(gca,'xtick',[],'ytick',1,'yticklabel','stress')
colormap('gray');freezeColors(gca);

axes('position',[0.1,0.96,0.85,0.01])
imagesc(cent_norm(log2(sorted_total_mol_1)))
set(gca,'xtick',[],'ytick',1,'yticklabel','tot. mol.')
colormap('jet');freezeColors(gca);

eval(['export_fig onlyneurons_hypoth_backspinv2_lev',num2str(splitlev),'_',date,'.pdf'])


% % % % % % % % % % % % % % % % % % % % % % % % % % %
T_cells_tmp_uni = unique(T_cells_tmp);
pmat_tmp = zeros(length(dataout_sorted_all_n(:,1)),length(T_cells_tmp_uni));
gr_center = zeros(length(T_cells_tmp_uni),1);
% for jjj=1:length(T_cells_tmp_uni)
%     jjj
%     gr_center(jjj) = mean(find(T_cells_tmp==T_cells_tmp_uni(jjj)));
%     [~,p] = ttest2(log2(dataout_sorted_all_n(:,T_cells_tmp==T_cells_tmp_uni(jjj))+1)',log2(dataout_sorted_all_n(:,T_cells_tmp~=T_cells_tmp_uni(jjj))+1)','tail','right');
%     pmat_tmp(:,jjj) = p;
% end
% for jjj=1:length(T_cells_tmp_uni)
%     jjj
%     gr_center(jjj) = mean(find(T_cells_tmp==T_cells_tmp_uni(jjj)));
%     p1 = zeros(length(dataout_sorted_all_n(:,1)),1);
%     inanal = mean(dataout_sorted_all_n(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2)>1;
%     for jjjj = 1:length(T_cells_tmp_uni)
%         if jjjj~=jjj
%             p = zeros(size(p1));
%             [~,p(inanal)] = ttest2(log2(dataout_sorted_all_n(inanal,T_cells_tmp==T_cells_tmp_uni(jjj))+1)'...
%                 ,log2(dataout_sorted_all_n(inanal,T_cells_tmp==T_cells_tmp_uni(jjjj))+1)','tail','right');            
%             p1(p>p1) = p(p>p1);
%         end
%     end
%     p1(~inanal) = 1;
%     pmat_tmp(:,jjj) = p1;
% end
sumpergene = sum(dataout_sorted_all_n,2);
molenrich_mat = zeros(length(dataout_sorted_all_n(:,1)),length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    gr_center(jjj) = mean(find(T_cells_tmp==T_cells_tmp_uni(jjj)));
    molenrich_mat(:,jjj) = sum(dataout_sorted_all_n(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2)./sumpergene*length(T_cells_tmp)/sum(T_cells_tmp==T_cells_tmp_uni(jjj));
end
[~,grord] = sort(gr_center);
molenrich_mat = molenrich_mat(:,grord);
[~,xi] = sort(molenrich_mat,'descend');


% [~,xi] = sort(pmat_tmp);
% [~,grord] = sort(gr_center);
% pmat_tmp = pmat_tmp(:,grord);
% [~,xi] = sort(pmat_tmp);

stress_ratio = zeros(length(T_cells_tmp_uni),1);
stress_ratio_p = zeros(length(T_cells_tmp_uni),1);
gr_stress_pval = ones(length(geneuni_1),length(T_cells_tmp_uni));
gr_stress_pval_n = ones(length(neurons_genes),length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    stress_ratio(jjj) = sum(sorted_stress_1(find(T_cells_tmp==T_cells_tmp_uni(jjj)))==1)/sum(T_cells_tmp==T_cells_tmp_uni(jjj));
    stress_ratio_p(jjj) = 1-binocdf( sum(sorted_stress_1(find(T_cells_tmp==T_cells_tmp_uni(jjj)))==1) -1 , sum(T_cells_tmp==T_cells_tmp_uni(jjj)), sum(sorted_stress_1==1)/length(sorted_stress_1)) ;
    if sum(T_cells_tmp==T_cells_tmp_uni(jjj) & sorted_stress_1==1)>1 & sum(T_cells_tmp==T_cells_tmp_uni(jjj) & sorted_stress_1==0)>1
        [~,p] = ttest2( log2(sorted_moldata_neurons_allgenes(:,T_cells_tmp==T_cells_tmp_uni(jjj) & sorted_stress_1==1)+1)',...
            log2(sorted_moldata_neurons_allgenes(:,T_cells_tmp==T_cells_tmp_uni(jjj) & sorted_stress_1==0)+1)','tail','right');
        gr_stress_pval(:,jjj) = p;
        
        [~,p] = ttest2( log2(dataout_sorted_all_n(:,T_cells_tmp==T_cells_tmp_uni(jjj) & sorted_stress_1==1)+1)',...
            log2(dataout_sorted_all_n(:,T_cells_tmp==T_cells_tmp_uni(jjj) & sorted_stress_1==0)+1)','tail','right');
        gr_stress_pval_n(:,jjj) = p;
    end
end
[~,xistress] = sort(gr_stress_pval);
[~,xistress_n] = sort(gr_stress_pval_n);
stress_genes_gr = cell(length(T_cells_tmp_uni),1);
for jjj=1:length(T_cells_tmp_uni)
    z=length(fdr_proc(gr_stress_pval(:,jjj),0.1));
    if z>3;
        disp(['gr ',num2str(T_cells_tmp_uni(jjj)),', ',num2str(z),...
            ' nonstress=',num2str(sum(T_cells_tmp==T_cells_tmp_uni(jjj) & sorted_stress_1==0)),', stress=',num2str(sum(T_cells_tmp==T_cells_tmp_uni(jjj) & sorted_stress_1==1))]);
        stress_genes_gr{jjj} = geneuni_1(fdr_proc(gr_stress_pval(:,jjj),0.1));
    end
end
stress_genes_vec = [];
for i=1:length(T_cells_tmp_uni)
    if ~isempty(stress_genes_gr{i})
        stress_genes_vec = [stress_genes_vec;stress_genes_gr{i}];
    end
end
stress_genes_vec_sort = sort(stress_genes_vec);
[stress_genes_vec_uni,ia] = unique(stress_genes_vec_sort,'first');
[~,ib] = unique(stress_genes_vec_sort,'last');



ind_gr_tmp_mark = xi(1,:);
ind_gr_tmp_mark = flipud(ind_gr_tmp_mark(:));
rmv = false(size(ind_gr_tmp_mark));
for i=2:length(ind_gr_tmp_mark)
    if ismember(ind_gr_tmp_mark(i),ind_gr_tmp_mark(1:i-1))
        rmv(i) = true;
    end
end
ind_gr_tmp_mark(rmv) = [];
gr_tmp_mark = neurons_genes(ind_gr_tmp_mark);
gr_tmp_mark = (gr_tmp_mark(:));

figure;
set(gcf,'position',[100,100,600,1000],'color','w')
barhight = 0.97/length(ind_gr_tmp_mark);
for jj=1:length(ind_gr_tmp_mark)
    disp([num2str(jj),',']);
    maxval = prctile(dataout_sorted_all_n(ind_gr_tmp_mark(jj),:),99.99);%0.9*max(data_tmpgr(ind_gr_tmp_mark(jj),:));
    n = length(T_cells_tmp);
    axes('position',[0.08,0.01+barhight*(jj-1),0.9,barhight])
    if mod(jj,2)==1;
        h = bar(dataout_sorted_all_n(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',1);hold on;axis tight;axis off;
    else
        h = bar(dataout_sorted_all_n(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',1,'facecolor',[200,0,0]/256);hold on;axis tight;axis off;
    end
    for jjj=1:length(cells_bor_tmpgr{lev})
        plot(cells_bor_tmpgr{lev}(jjj)*[1,1]-1,[0,maxval],'--','linewidth',0.5,'color',[200,200,200]/256)
    end
    set(gca,'xlim',[0,n],'ylim',[0,maxval])
    text(-0.5,maxval/2,gr_tmp_mark{jj},'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',7)
end

eval(['export_fig only_neurons_2nd_barplots_molenrich_sum_backspinv2_lev',num2str(splitlev),'_',date,'.pdf'])

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
colorvec2 = rand(max(T_cells_tmp),3);
tmp_Tcells = zeros(length(T_cells_tmp),1,3);
for i=1:max(T_cells_tmp)
   tmp_Tcells(T_cells_tmp==i,1,:) = repmat(colorvec2(i,:),sum(T_cells_tmp==i),1);
end
tmp_Tcells = reshape(tmp_Tcells,1,length(tmp_Tcells),3);


mol_data_sorted = dataout_sorted_n;
ordtmp = 'ABCD';
foldername = 'backspinv2_neurons_lev7_june4_2015/';
figure('color','w','position',[200,200,650,1000],'visible','on');
k=0;


for i=1:length(mol_data_sorted(:,1))
    i
    axes('position',[0.1,0.97-(mod(i-1,10)+1)*0.091,0.85,0.091])
    tmp = mol_data_sorted(i,:);
    bar(tmp);hold on;
    plot([find(diff(T_cells_tmp')~=0)*[1,1]]',[repmat([0,max(tmp)],sum(diff(T_cells_tmp)~=0),1)]','--','color',get_RGB('grey'))
    axis tight
    text(10,max(tmp),[neurons_genes_sorted{i},', i=',num2str(i)],'VerticalAlignment','top');
    if mod(i-1,10)+1==10
        set(gca,'xtick',[0:100:length(T_cells_tmp(1,:))],'fontsize',6);%,'position',[0.1,0.98-(mod(i,10)+1)*0.092,0.85,0.092]);
        axes('position',[0.1,0.01,0.85,0.0075])
        imagesc(sorted_age_1)
        axis off
        axes('position',[0.1,0.0175,0.85,0.0075])
        imagesc(cent_norm(log2(sorted_total_mol_1)))
        axis off
        axes('position',[0.1,0.025,0.85,0.0075])
        imagesc(cent_norm(sorted_stress_1))
        axis off        
        axes('position',[0.1,0.0325,0.85,0.0075])
        image(tmp_Tcells); hold on;
        plot([find(diff(T_cells_tmp')~=0)*[1,1]]',[repmat([0,max(tmp)],sum(diff(T_cells_tmp)~=0),1)]','k')        
        axis off
        
        
        axes('position',[0.1,0.97,0.85,0.0075])
        image(tmp_Tcells); hold on;
        plot([find(diff(T_cells_tmp')~=0)*[1,1]]',[repmat([0,max(tmp)],sum(diff(T_cells_tmp)~=0),1)]','k')        
        axis off                
        axes('position',[0.1,0.9775,0.85,0.0075])
        imagesc(cent_norm(log2(sorted_total_mol_1)))
        axis off
        axes('position',[0.1,0.985,0.85,0.0075])
        imagesc(cent_norm(sorted_stress_1))
        axis off          
        axes('position',[0.1,0.9925,0.85,0.0075])
        imagesc(sorted_age_1)
        axis off    
        
        k=k+1;
        fliplr([mod(k/1,4),mod(k/4,4),mod(k/4^2,4),mod(k/4^3,4),mod(k/4^4,4),mod(k/4^5,4)]);
%         save2pdf([foldername,ordtmp(floor(ans)+1),'_barplot_gene_hypoth_heatmaporder_',num2str(k),'.pdf'],gcf,300)
        eval(['export_fig ',foldername,ordtmp(floor(ans)+1),'_barplot_gene_hypoth_heatmaporder_',num2str(k),'.pdf'])
        close all
        figure('color','w','position',[200,200,650,1000],'visible','off');
    else
        set(gca,'xtick',[10:100:length(T_cells_tmp(1,:))],'xticklabel',cell(length([10:50:length(T_cells_tmp(1,:))]),1),'fontsize',8);%...
%             ,'position',[0.1,0.98-(mod(i,10)+1)*0.092,0.85,0.092]);
    end
end
if mod(i-1,10)+1~=10
    set(gca,'xtick',[0:100:length(T_cells_tmp(1,:))],'fontsize',6);%,'position',[0.1,0.98-(mod(i,10)+1)*0.092,0.85,0.092]);
        axes('position',[0.1,0.01,0.85,0.0075])
        imagesc(sorted_age_1)
        axis off
        axes('position',[0.1,0.0175,0.85,0.0075])
        imagesc(cent_norm(log2(sorted_total_mol_1)))
        axis off
%         axes('position',[0.1,0.025,0.85,0.0075])
%         imagesc(cent_norm(sorted_strees6h_1))
%         axis off        
        axes('position',[0.1,0.0325,0.85,0.0075])
        image(tmp_Tcells); hold on;
        plot([find(diff(T_cells_tmp')~=0)*[1,1]]',[repmat([0,max(tmp)],sum(diff(T_cells_tmp)~=0),1)]','k')        
        axis off
        
        
        axes('position',[0.1,0.97,0.85,0.0075])
        image(tmp_Tcells); hold on;
        plot([find(diff(T_cells_tmp')~=0)*[1,1]]',[repmat([0,max(tmp)],sum(diff(T_cells_tmp)~=0),1)]','k')        
        axis off                
        axes('position',[0.1,0.9775,0.85,0.0075])
        imagesc(cent_norm(log2(sorted_total_mol_1)))
        axis off
%         axes('position',[0.1,0.985,0.85,0.0075])
%         imagesc(cent_norm(sorted_strees6h_1))
%         axis off          
        axes('position',[0.1,0.9925,0.85,0.0075])
        imagesc(sorted_age_1)
        axis off        
    k=k+1;
    fliplr([mod(k/1,4),mod(k/4,4),mod(k/4^2,4),mod(k/4^3,4),mod(k/4^4,4),mod(k/4^5,4)]);
%     save2pdf([foldername,ordtmp(floor(ans)+1),'_barplot_gene_hypoth_heatmaporder_',num2str(k),'.pdf'],gcf,300)
    eval(['export_fig ',foldername,ordtmp(floor(ans)+1),'_barplot_gene_hypoth_heatmaporder_',num2str(k),'.pdf'])
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



sumpergene = sum(dataout_sorted_all_n,2);
meanpergene = mean(dataout_sorted_all_n,2);%
molenrich_mat = zeros(length(dataout_sorted_all_n(:,1)),length(T_cells_tmp_uni));
meangr_mat = zeros(length(dataout_sorted_all_n(:,1)),length(T_cells_tmp_uni));
meangrpos_mat = zeros(length(dataout_sorted_all_n(:,1)),length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    gr_center(jjj) = mean(find(T_cells_tmp==T_cells_tmp_uni(jjj)));
    molenrich_mat(:,jjj) = mean(dataout_sorted_all_n(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2)./meanpergene;
    meangrpos_mat(:,jjj) = mean(dataout_sorted_all_n(:,T_cells_tmp==T_cells_tmp_uni(jjj))>0,2);
end
[~,grord] = sort(gr_center);
meangrpos_mat = meangrpos_mat(:,grord);
molenrich_mat = molenrich_mat(:,grord);
[~,xi0] = sort(molenrich_mat.*meangrpos_mat.^0,'descend');
[~,xi0p5] = sort(molenrich_mat.*meangrpos_mat.^1,'descend');
[~,xi1] = sort(molenrich_mat.*meangrpos_mat.^2,'descend');


ind_gr_tmp_mark = [xi0(1:5,:);xi0p5(1:5,:);xi1(1:5,:)];
ind_gr_tmp_mark = flipud(ind_gr_tmp_mark(:));
rmv = false(size(ind_gr_tmp_mark));
for i=2:length(ind_gr_tmp_mark)
    if ismember(ind_gr_tmp_mark(i),ind_gr_tmp_mark(1:i-1))
        rmv(i) = true;
    end
end
ind_gr_tmp_mark(rmv) = [];
gr_tmp_mark = neurons_genes(ind_gr_tmp_mark);
gr_tmp_mark = (gr_tmp_mark(:));

datamarkers = dataout_sorted_all_n(ind_gr_tmp_mark,:);
datamarkers_cn = cent_norm(log2(datamarkers+1));


figure;
set(gcf,'position',[100,100,600,1000],'color','w')
axes('position',[0.1,0.02,0.88,0.95])
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
hold on;
linewid =0.5;
bor_color = 'grey11';%'green1';%
for jj=1:length(cells_bor_tmpgr{lev})
    plot(cells_bor_tmpgr{lev}(jj)*[1,1]-0.5,[1,length(gr_tmp_mark)],'-','linewidth',linewid,'color',get_RGB(bor_color))
end
set(gca,'xtick',[0:1000:length(cell_neurons_sorted)],'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark)
colormap('summer');
freezeColors(gca);


table = m2c([T_cells_tmp;datamarkers]');
table = [[{'','cluster'},gr_tmp_mark'];[cell_neurons_sorted,table]];

saveCellFile(table,['markertable_hypo_neuronsonly_backspinv2_lev',num2str(splitlev),'_',date,'.txt'])











