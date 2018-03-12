tic
clear all
close all
moldata_all = [];
source_plate = [];
cellid_all = [];
load /data2/WG9600_stuff/first_data_test/WG960016_NOsingletons_ctx_test_24-Jun-2016.mat
moldata_all = [moldata_all,moldata];
cellid_all = [cellid_all,cellid];
source_plate = [source_plate,repmat({'960016'},1,length(moldata(1,:)))];
load /data2/WG9600_stuff/first_data_test/WG960015_NOsingletons_ctx_test_24-Jun-2016.mat
moldata_all = [moldata_all,moldata];
cellid_all = [cellid_all,cellid];
source_plate = [source_plate,repmat({'960015'},1,length(moldata(1,:)))];
load /data2/WG9600_stuff/first_data_test/WG960017_NOsingletons_ctx_test_27-Jun-2016.mat
moldata_all = [moldata_all,moldata];
cellid_all = [cellid_all,cellid];
source_plate = [source_plate,repmat({'960017'},1,length(moldata(1,:)))];
load /data2/WG9600_stuff/first_data_test/WG960011_ctx_test_27-Jun-2016.mat
moldata_all = [moldata_all,moldata];
cellid_all = [cellid_all,cellid];
source_plate = [source_plate,repmat({'960011'},1,length(moldata(1,:)))];
load /data2/WG9600_stuff/first_data_test/WG960021_ctx_2lanes_test_19-Jul-2016.mat
moldata_all = [moldata_all,moldata];
cellid_all = [cellid_all,cellid];
source_plate = [source_plate,repmat({'960021'},1,length(moldata(1,:)))];


moldata = moldata_all;
cellid = cellid_all;

a = loadCellFile_turbo('/data2/WG9600_stuff/first_data_test/markertable_WG9600_11_15_16_17_21_ctx_lev7_01-Sep-2016.txt',1);
cellid_class_lev1 = a(2:end,[1,4]);
cellid_class_lev1( cellfun(@isempty,cellid_class_lev1(:,2)),:) = [];

cellid_neurons = cellid_class_lev1(strcmpi(cellid_class_lev1(:,2),'neurons'),1);
cellid_nonneurons = cellid_class_lev1(~strcmpi(cellid_class_lev1(:,2),'neurons'),1);
[~,locneurons] = ismember(cellid_neurons,cellid);
[~,locnonneurons] = ismember(cellid_nonneurons,cellid);

[~,ptt] = ttest2(log2(moldata(:,locneurons)+1)', log2(moldata(:,locnonneurons)+1)', 'tail', 'right');
neurons_genes_ind = fdr_proc(ptt,0.01);

moldata = moldata(neurons_genes_ind,locneurons);
cellid = cellid(locneurons);
source_plate = source_plate(locneurons);
genes_sym = genes_sym(neurons_genes_ind);

tic
splitlev = 6;
Nfeture = 200;
Nfeture1 = 800;
N_to_backspin = 10;
N_to_cells = 500;
mean_tresh = 0.01;
fdr_th = 0.3;
min_gr_cells = 5;
min_gr_genes = 20;
stop_th = [0.5,0.5];
flag_val_stip = 1;
% permute the cells order
rng(1);
cell_perm = randperm(length(cellid));
gene_perm = randperm(length(genes_sym));
genes_sym = genes_sym(gene_perm);
cellid = cellid(cell_perm);
source_plate = source_plate(cell_perm);
moldata = moldata(gene_perm,cell_perm);

[dataout_sorted,genes_order,cells_order,genes_gr_level,cells_gr_level,genes_bor_level,cells_bor_level] =...
    backSpinSplit_v7(moldata,splitlev,Nfeture,Nfeture1,N_to_backspin,N_to_cells,mean_tresh,fdr_th,min_gr_cells,min_gr_genes,stop_th,flag_val_stip);
toc




data_sorted_all = moldata(:,cells_order);
sorted_data_cn = cent_norm(log2_center(dataout_sorted));
cellid_sorted = cellid(cells_order);
source_plate_sorted = source_plate(cells_order);


T_cells_tmp = cells_gr_level(:,splitlev+1)';
T_genes_tmp = genes_gr_level(:,splitlev+1)';

cells_bor_tmpgr = cells_bor_level;
genes_bor_tmpgr = genes_bor_level;

lev = splitlev;
cells_bor_tmpgr{lev} = [cells_bor_tmpgr{lev};length(T_cells_tmp)+1];
genes_bor_tmpgr{lev} = [genes_bor_tmpgr{lev};length(T_genes_tmp)+1];

figure;
set(gcf,'position',[100,100,450,750],'color','w')
axes('position',[0.1,0.02,0.88,0.85])
imagesc(sorted_data_cn,[prctile(sorted_data_cn(:),1),prctile(sorted_data_cn(:),99)]);
hold on;
linewid =0.5;
bor_color = 'grey11';%'green1';%
for jj=1:length(cells_bor_tmpgr{lev})
    if jj>1
        plot(cells_bor_tmpgr{lev}(jj)*[1,1]-0.5,[genes_bor_tmpgr{lev}(jj-1),genes_bor_tmpgr{lev}(jj)],'--','linewidth',linewid,'color',get_RGB(bor_color))
        plot([cells_bor_tmpgr{lev}(jj-1),cells_bor_tmpgr{lev}(jj)]-0.5,genes_bor_tmpgr{lev}(jj)*[1,1] ,'--','linewidth',linewid,'color',get_RGB(bor_color))
        plot(cells_bor_tmpgr{lev}(jj-1)*[1,1]-0.5,[genes_bor_tmpgr{lev}(jj-1),genes_bor_tmpgr{lev}(jj)]-0.5,'--','linewidth',linewid,'color',get_RGB(bor_color))
        plot([cells_bor_tmpgr{lev}(jj-1)-0.5,cells_bor_tmpgr{lev}(jj)],genes_bor_tmpgr{lev}(jj-1)*[1,1]-0.5 ,'--','linewidth',linewid,'color',get_RGB(bor_color))
    else
        plot(cells_bor_tmpgr{lev}(jj)*[1,1]-0.5,[1,genes_bor_tmpgr{lev}(jj)]-0.5,'--','linewidth',linewid,'color',get_RGB(bor_color))
        plot([1,cells_bor_tmpgr{lev}(jj)]-0.5,genes_bor_tmpgr{lev}(jj)*[1,1] -0.5,'--','linewidth',linewid,'color',get_RGB(bor_color))
        plot(1*[1,1]-0.5,[1,genes_bor_tmpgr{lev}(jj)]-0.5,'--','linewidth',linewid,'color',get_RGB(bor_color))
        plot([1,cells_bor_tmpgr{lev}(jj)]-0.5,1*[1,1] -0.5,'--','linewidth',linewid,'color',get_RGB(bor_color))
    end
end
set(gca,'xtick',[0:200:length(cells_order)],'ytick',[0:100:length(genes_order)], 'fontsize', 6)
colormap('summer');
freezeColors(gca);
%eval(['export_fig backspin_DI_WG96007_cortex_',date,'.pdf'])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

T_cells_tmp_uni = unique(T_cells_tmp);
gr_center = zeros(length(T_cells_tmp_uni),1);

% % % % % % % % % % % % % % %
% Marker genes barplots


meanpergene = mean(data_sorted_all,2);%mean(log2(data_sorted_all+1),2);
molenrich_mat = zeros(length(data_sorted_all(:,1)),length(T_cells_tmp_uni));
meangr_mat = zeros(length(data_sorted_all(:,1)),length(T_cells_tmp_uni));
meangrpos_mat = zeros(length(data_sorted_all(:,1)),length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    gr_center(jjj) = mean(find(T_cells_tmp==T_cells_tmp_uni(jjj)));
%     molenrich_mat(:,jjj) = sum(data_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2)./sumpergene*length(T_cells_tmp)/sum(T_cells_tmp==T_cells_tmp_uni(jjj));
    molenrich_mat(:,jjj) = mean(data_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2)./meanpergene;
%     meangr_mat(:,jjj) = mean(log2(data_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj))+1),2);
    meangrpos_mat(:,jjj) = mean(data_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj))>0,2);
end
% molenrich_mat(meangrpos_mat<0.1) = 0;
meangrpos_mat(meangrpos_mat<0.2) = 0;
molenrich_mat(meanpergene==0 | isnan(meanpergene),:) = 0;
molenrich_mat = molenrich_mat.*meangrpos_mat.^1;
[~,grord] = sort(gr_center);
molenrich_mat = molenrich_mat(:,grord);
[~,xi] = sort(molenrich_mat,'descend');


ind_gr_tmp_mark = xi(1:3,:);
ind_gr_tmp_mark = flipud(ind_gr_tmp_mark(:));
rmv = false(size(ind_gr_tmp_mark));
for i=2:length(ind_gr_tmp_mark)
    if ismember(ind_gr_tmp_mark(i),ind_gr_tmp_mark(1:i-1))
        rmv(i) = true;
    end
end
ind_gr_tmp_mark(rmv) = [];
gr_tmp_mark = genes_sym(ind_gr_tmp_mark);
gr_tmp_mark = (gr_tmp_mark(:));

figure;
set(gcf,'position',[100,100,450,1000],'color','w')


barhight = 0.96/length(ind_gr_tmp_mark);
for jj=1:length(ind_gr_tmp_mark)
    disp([num2str(jj),','])
    maxval = prctile(data_sorted_all(ind_gr_tmp_mark(jj),:),99.99);%0.9*max(data_tmpgr(ind_gr_tmp_mark(jj),:));
    n = length(T_cells_tmp);
    axes('position',[0.1,0.01+barhight*(jj-1),0.85,barhight])
    %     if mod(jj,2)==1;
    %         bar(n/2,maxval,'facecolor',get_RGB('grey81'),'edgecolor','none','BarWidth',n);hold on;axis tight;axis off;
    %     end
    if mod(jj,2)==1;
        h = bar(data_sorted_all(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',2);hold on;axis tight;axis off;
    else
        h = bar(data_sorted_all(ind_gr_tmp_mark(jj),:),'edgecolor','none','BarWidth',2,'facecolor',[200,0,0]/256);hold on;axis tight;axis off;
    end
    for jjj=1:length(cells_bor_tmpgr{lev})
        plot(cells_bor_tmpgr{lev}(jjj)*[1,1]-1,[0,maxval],'--','linewidth',0.5,'color',[200,200,200]/256); hold on;
    end
    
    set(gca,'xlim',[0,n],'ylim',[0,maxval], 'fontsize', 6)
    text(-0.5,maxval/2,gr_tmp_mark{jj},'HorizontalAlignment','right','VerticalAlignment','middle','fontsize',7)
end
%eval(['export_fig marker_barplots_DI_WG96007_cortex_',date,'.pdf'])

eval(['export_fig barplots_neuronsonly_WG960011_15_16_17_21_lev',num2str(splitlev),'_',date,'.pdf'])

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



sumpergene = sum(data_sorted_all,2);
meanpergene = mean(data_sorted_all,2);%
molenrich_mat = zeros(length(data_sorted_all(:,1)),length(T_cells_tmp_uni));
meangr_mat = zeros(length(data_sorted_all(:,1)),length(T_cells_tmp_uni));
meangrpos_mat = zeros(length(data_sorted_all(:,1)),length(T_cells_tmp_uni));
for jjj=1:length(T_cells_tmp_uni)
    jjj
    gr_center(jjj) = mean(find(T_cells_tmp==T_cells_tmp_uni(jjj)));
    molenrich_mat(:,jjj) = mean(data_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj)),2)./meanpergene;
    meangrpos_mat(:,jjj) = mean(data_sorted_all(:,T_cells_tmp==T_cells_tmp_uni(jjj))>0,2);
end
[~,grord] = sort(gr_center);
meangrpos_mat(meangrpos_mat<0.2) = 0;
molenrich_mat(meanpergene==0 | isnan(meanpergene),:) = 0;
meangrpos_mat = meangrpos_mat(:,grord);
molenrich_mat = molenrich_mat(:,grord);
[~,xi0] = sort(molenrich_mat.*meangrpos_mat.^0.001,'descend');
[~,xi0p5] = sort(molenrich_mat.*meangrpos_mat.^0.5,'descend');
[~,xi1] = sort(molenrich_mat.*meangrpos_mat.^1,'descend');


ind_gr_tmp_mark = [xi0(1:15,:);xi0p5(1:15,:);xi1(1:15,:)];
ind_gr_tmp_mark = flipud(ind_gr_tmp_mark(:));
rmv = false(size(ind_gr_tmp_mark));
for i=2:length(ind_gr_tmp_mark)
    if ismember(ind_gr_tmp_mark(i),ind_gr_tmp_mark(1:i-1))
        rmv(i) = true;
    end
end
ind_gr_tmp_mark(rmv) = [];
gr_tmp_mark = genes_sym(ind_gr_tmp_mark);
gr_tmp_mark = (gr_tmp_mark(:));
molenrich_mat_mark = zeros(length(ind_gr_tmp_mark),length(T_cells_tmp_uni));
meangrpos_mat_mark = zeros(length(ind_gr_tmp_mark),length(T_cells_tmp_uni));
order_gr = [T_cells_tmp(diff(T_cells_tmp)~=0),T_cells_tmp(end)];
for jjj=1:length(T_cells_tmp_uni)
    jjj
    molenrich_mat_mark(:,jjj) = mean(data_sorted_all(ind_gr_tmp_mark,T_cells_tmp==order_gr(jjj)),2)./meanpergene(ind_gr_tmp_mark);
    meangrpos_mat_mark(:,jjj) = mean(data_sorted_all(ind_gr_tmp_mark,T_cells_tmp==order_gr(jjj))>0,2);
end
[~,imax] = max(molenrich_mat_mark.*meangrpos_mat_mark.^0.5,[],2);
[~,xi] = sort(imax);
ind_gr_tmp_mark = ind_gr_tmp_mark(xi);

datamarkers = data_sorted_all(ind_gr_tmp_mark,:);
datamarkers_cn = cent_norm(log2(datamarkers+1));


figure;
set(gcf,'position',[100,100,450,1000],'color','w')
axes('position',[0.1,0.02,0.88,0.96])
imagesc(datamarkers_cn,[prctile(datamarkers_cn(:),1),prctile(datamarkers_cn(:),99)]);
hold on;
linewid =0.5;
bor_color = 'grey11';%'green1';%
for jj=1:length(cells_bor_tmpgr{lev})
    plot(cells_bor_tmpgr{lev}(jj)*[1,1]-0.5,[1,length(gr_tmp_mark)],'-','linewidth',linewid,'color',get_RGB(bor_color))
end
set(gca,'xtick',[0:200:length(cells_order)],'ytick',[1:length(gr_tmp_mark)],'yticklabel',gr_tmp_mark(xi), 'fontsize', 8)
colormap('summer');
freezeColors(gca);

eval(['export_fig markertable_neuronsonly_WG960011_15_16_17_21_lev',num2str(splitlev),'_',date,'.pdf'])


table1 = m2c([T_cells_tmp;datamarkers]');
table1 = [[{'','tissue','cluster'},gr_tmp_mark(xi)'];[cellid_sorted',source_plate_sorted',table1]];

saveCellFile(table1,['markertable_neurons_WG9600_11_15_16_17_21_ctx_lev',num2str(splitlev),'_',date,'.txt'])



% % % % % % % % % % % % % % % % % % % % % % %
% 
% m_v = mean(moldata,2);
% cv_v = std(moldata,[],2)./m_v;
% log2_m = log2(m_v);
% log2_cv = log2(cv_v);
% x0 = [-0.5,1];
% [param_fit,fval] =  run_min_cvfit(log2_m,log2_cv,x0);
% param_fit = round(param_fit*100)/100;
% log2_cv_fit = log2( (2.^log2_m).^param_fit(1) + param_fit(2));
% tmp = log2_cv - log2_cv_fit;
% [~,xi1] = sort(tmp,'descend');
% corr_filt = xi1(1:1000);
% 
% moldata_norm = moldata./ repmat(sum(moldata),length(moldata(:,1)),1);
% data = log2_center(moldata_norm(corr_filt,:));
% % data = log2_center(moldata(corr_filt,:));
% 
% indgood = 1:length(data(1,:));%find(T_cells_tmp>0); %
% 
% no_dims = 2;
% initial_dims = 20;%length(gr_tmp_mark);%20;
% perplexity = 20;
% epsilon = 100;
% dist_flag = 2;
% max_iter = 1000;
% 
% tic
% mappedX_cell = tsne((data(:,indgood))', [], no_dims, initial_dims, perplexity,epsilon,dist_flag,max_iter);
% toc
% % % % %
% 
% figure;
% set(gcf,'position',[100,100,900,900],'color','w')
% plot(mappedX_cell(:,1),mappedX_cell(:,2),'.');
% 
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% 
% list = {'Gad1','Npy','Rasgrf2','Pvrl3','Rorb','Rprm','Foxp2','Cplx3','Syt6'};
% figure;
% set(gcf,'color','w','position',[20,20,1000,1500])
% for i=1:length(list)
%     genePlot = list{i};
%     markergene = (moldata(strcmpi(genes_sym,genePlot),:));
%     if max(markergene)>prctile(markergene,90) & prctile(markergene,90)~=0
%         markergene(markergene>prctile(markergene,99)) = prctile(markergene,99);
%         markergene(markergene<prctile(markergene,3)) = prctile(markergene,3);
%     end
%     c_rgb = [1,0,0];rand([1,3]);
%     %     markergene_color = [interp1([min(markergene),max(markergene)],[0,1],markergene'),zeros(size(markergene'))...
%     %         ,interp1([min(markergene),max(markergene)],[1,0],markergene')];
%     markergene_color = [interp1([min(markergene),max(markergene)],[0.7,c_rgb(1)],markergene'),...
%         interp1([min(markergene),max(markergene)],[0.7,c_rgb(2)],markergene')...
%         ,interp1([min(markergene),max(markergene)],[0.7,c_rgb(3)],markergene')];
%     subplot(3,3,i);
%     scatter(mappedX_cell(:,1),mappedX_cell(:,2),80,markergene_color,'.')
%     set(gca,'xlim',[-150,150],'ylim',[-150,150])
%     title(genePlot);
%     axis tight
% end


