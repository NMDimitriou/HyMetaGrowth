function plot_hist_z(coord,ca_coord1,tp,name)

gc=struct2array(coord);
gc=struct2table(gc);
gc=table2array(gc);
gc_ca=struct2array(ca_coord1);
gc_ca=struct2table(gc_ca);
gc_ca=table2array(gc_ca);
[l1,d1]=size(gc);
[l2,d2]=size(gc_ca);
%h=15;

for i=1:length(tp)
    
    gc_z=[];
    gc_ca_z=[];
    
    for j=1:l1       
        gc_z = [gc_z; gc{j,i}];      
    end
%    lc = length(gc_z);
%    gc_z(:,3)=gc_z(:,3)+sqrt(h)*rand(lc,1);
    
    for j=1:l2
        gc_ca_z = [gc_ca_z; gc_ca{j,i}];
    end
    
    f1=figure('Position', [100, 100, 800, 700],'Visible','off');
    hold on
    histogram(gc_z(:,3),100,'FaceAlpha',0.5,'DisplayName','Experiments','Normalization','probability')  
    histogram(gc_ca_z(:,3),100,'FaceAlpha',0.5,'DisplayName','Simulations','Normalization','probability')
    legend
    title(tp(i))
    xlabel('z ($\mu m$)')
    axis([0 900 0 inf])
    hold off
    set(f1,'Units','Inches');
    pos = get(f1,'Position');
    set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    disp('Saving...')
    print(f1,[name tp{i} '.png'],'-r350','-dpng')
    
end

end