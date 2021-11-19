function [samp1,sampc1,coord,count,ca_coord1] = import_coord(nsamp,tp,ntp,file1,file2,group1,gname1,groupcs1,gname_ca1,...
          samp1,sampc1,coord,count,ca_coord1,k)

for i=1:nsamp
   
    for j=1:ntp
       
        samp1{i+k,j}=dir([file1 tp{j} '/' group1{i}]); %
        coord.(gname1{i}).(tp{j})=readmatrix([samp1{i+k,j}.folder '/' samp1{i+k,j}.name]);
        count.(gname1{i}).(tp{j})=length(coord.(gname1{i}).(tp{j})(:,1));
        sampc1{i+k,j}=dir([file2 tp{j} '/' groupcs1{i}]);
        ca_coord1.(gname_ca1{i}).(tp{j})=readmatrix([sampc1{i+k,j}.folder '/' sampc1{i+k,j}.name]);
    end
    %sampc1{i}=dir([file2 groupc1{i}]);
    %ca_coord1.(gname_ca1{i}) = readmatrix([sampc1{i}.folder '/'  sampc1{i}.name]);
end

end