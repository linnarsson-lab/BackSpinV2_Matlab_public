function labels = knn_corr_classi(data,class_key,numneibor)

class_key_uni = unique(class_key);
if length(data(:,1))>=2
    if mod(numneibor,2)==0
        numneibor = numneibor+1;
    end
    
    D = corr_mat(data);
    [~, ind] = sort(D,'descend');
    kneighbours = ind(2:1+numneibor, :);
    n_class = class_key(kneighbours);
    if length(n_class(1,:))~=length(data(1,:))
        n_class = [n_class(:)]';
    end
    vote_perclass = zeros(length(class_key_uni),length(data(1,:)));
    for i=1:length(class_key_uni)
        vote_perclass(i,:) = sum(n_class==class_key_uni(i),1);
    end
    
    [~,indmax] = max(vote_perclass);
    labels = class_key_uni(indmax);
else
    labels = class_key;
end






