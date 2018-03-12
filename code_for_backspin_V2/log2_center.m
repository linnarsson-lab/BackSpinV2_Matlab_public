function datalog = log2_center(data)
datalog = log2(data+1);
datalog = datalog - repmat(mean(datalog,2),1,length(data(1,:)));