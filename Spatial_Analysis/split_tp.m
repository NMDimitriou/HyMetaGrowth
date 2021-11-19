function ca_coord1_tsp = split_tp(groupc1, gname1, ca_coord1, ts)

for i=1:length(groupc1)
    B = [];
    c=1;
    len=length(ca_coord1.(gname1{i})(:,1));
    for j=2:len

        if(ca_coord1.(gname1{i})(j,1)==ca_coord1.(gname1{i})(j-1,1))
                B=[B;ca_coord1.(gname1{i})(j-1,:)];
                ca_coord1_tsp.(gname1{i}).(ts{c})=B;
        else
                c=c+1;
                B=[];
                ca_coord1_tsp.(gname1{i}).(ts{c})=B;
        end
    end
end

end