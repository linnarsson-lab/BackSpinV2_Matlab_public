function D=corr_mat(A)


[n,m]=size(A);

c = cov(A);
D = corrcov(c,1);
% D=zeros(m);
% 
% for i=1:m
%     D(i,i) = 1;
%     for j=i+1:m
%         R = corrcoef(A(:,i),A(:,j));
%         D(i,j) = R(1,2);
% %         R = corr(A(:,i),A(:,j));
% %         D(i,j) = R(1,1);
%         D(j,i) = D(i,j);
%     end
% end











