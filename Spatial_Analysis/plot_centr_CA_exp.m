function plot_centr_CA_exp(groupc1,ca_coord1_tsp,coord1,gname1,tp,name)

for i=1:length(groupc1)
    
%    f=figure('Name',gname1{i},'Position', [100, 100, 1400, 1200]);
    NameArray = {'MarkerSize'};
    ValueArray = {3};
    for j=1:length(tp)
        A=ca_coord1_tsp.(gname1{i}).(tp{j});
        B=coord1.(gname1{i}).(tp{j});
        f=figure('Name',gname1{i},'Position', [100, 100, 800, 600]);
%        h(j)=subplot(3,3,j);
%        hax=gca;
%        hax.FontSize=16;
%        cmap = jet(length(A(:,3)));
         AxesH = axes('Parent', f, ...
             'NextPlot', 'add');  % Equivalent to: "hold on"
         xlabel('x $(\mu m)$','FontSize',30);
         ylabel('y $(\mu m)$','FontSize',30);
         zlabel('z $(\mu m)$','FontSize',30);
         title(tp{j},'Interpreter','latex','FontSize',35);
         axis([0 2500 0 2500 0 917])

%        for k=1:length(A(:,3))
%            plot3(A(k,1),A(k,2),A(k,3),'o','color',cmap(k,:),'markerfacecolor',cmap(k,:))
%        end
        scatter3(A(:,4),A(:,5),A(:,6),5,'filled','MarkerFaceColor','b');
        scatter3(B(:,1),B(:,2),B(:,3),5,'filled' ,'MarkerFaceColor','r');
        view(50,35)
        colormap(jet);
%        A(j)=copyobj(allchild(get(A(j),'CurrentAxes')),h(j));
%        view(20,20)
        %view(2)
%        set(A(j),NameArray,ValueArray)
      
        %hold off
        clear A
        clear B
        set(f,'Units','Inches');
        pos = get(f,'Position');
        set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
        disp('Saving...')
        print(f,[name '_' gname1{i} '_tp_' num2str(tp{j}) '.png'],'-r350','-dpng')
    end
    %subplot(3,3,8);
    %hcx=gca;
    %hcx.FontSize=16;
    %hold on
    %plot(time,struct2array(count.(gname{i})),'.-','LineWidth',4,'MarkerSize',30);
    %errorbar(time,gcm,gcstd,'.-','LineWidth',4,'MarkerSize',30)
    %xlabel('time $(day)$','Interpreter','latex')
    %ylabel('Nuclei count','Interpreter','latex')
    %hold off
    
%    set(f,'Units','Inches');
%    pos = get(f,'Position');
%    set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%    disp('Saving...')
%    print(f,['centroids_' gname1{i} '.png'],'-r350','-dpng')
end


end