function dy=d1calc(x,y)
[m,n]=size(x);
flip_flag=0;
if m>n
    x=x';
    y=y';
    n=m;
    flip_flag=1;
end

if length(unique(diff(x)))==1
    h=unique(diff(x));
    dy=[(y(2)-y(1))/h,(y(3:end)-y(1:end-2))/2/h,(y(end)-y(end-1))/h];
else
    dy=ones(1,n);
    dy(1)=(y(2)-y(1))/(x(2)-x(1));
    dy(end)=(y(end)-y(end-1))/(x(end)-x(end-1));
    dy(2:end-1)=(y(3:end)-y(1:end-2))./(x(3:end)-x(1:end-2));
end

if flip_flag==1
    dy=dy';
end



