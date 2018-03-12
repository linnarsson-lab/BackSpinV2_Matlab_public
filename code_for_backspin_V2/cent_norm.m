%input should be after log transform
function data_cent_norm=cent_norm(data);

data_mean=mean(data,2);
data_std=std(data')';
data_cent_norm=zeros(size(data));
for i=1:length(data(:,1))
    if data_std(i)>0
        data_cent_norm(i,:)=(data(i,:)-data_mean(i))/data_std(i);
    end
end