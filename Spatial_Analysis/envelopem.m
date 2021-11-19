function varargout = envelopem(A)
%% Plots the envelope of the summary function for a three-dimensional point pattern
% from all the samples.
% input: the matrix obtained from envelope.pp3 function of spatstat package
% for all samples
% of R
% output: the plot
% Author: Nikolaos M. Dimitriou, 
% McGill University, 2020

for k=1:length(A)
    
    Exist_obs = strcmp('obs',A{k}{1}.Properties.VariableNames);
    Exist_iso = strcmp('iso',A{k}{1}.Properties.VariableNames);
    
    if(Exist_obs(Exist_obs==1))
    
        x{k}=A{k}{1}.r';
        ytheo{k}=A{k}{1}.theo'-4*pi*x{k}.^3/3;

        for i=1:length(A{k})   
            yobs_all{k}(:,i) = A{k}{i}.obs;
            ylo_all{k}(:,i)  = A{k}{i}.lo;
            yhi_all{k}(:,i)  = A{k}{i}.hi;
        end

        yobs_mean{k} = mean(yobs_all{k},2)' -4*pi*x{k}.^3/3;
        yobs_std{k}  = std(yobs_all{k},0,2)'; %)./sqrt(length(yobs_all{k}(1,:))
        ylo_min{k}  = min(ylo_all{k},[],2)'  -4*pi*x{k}.^3/3;
        yhi_max{k}  = max(yhi_all{k},[],2)'  -4*pi*x{k}.^3/3;

        yobs_ms{k}   = [yobs_mean{k}-yobs_std{k}; yobs_mean{k}+yobs_std{k}];
        
    elseif(Exist_iso(Exist_iso==1))
        
        x{k}=A{k}{1}.r';
        ytheo{k}=A{k}{1}.iso'-4*pi*x{k}.^3/3;
        
        for i=1:length(A{k})   
            yobs_all{k}(:,i) = A{k}{i}.iso;
        end
        
        yobs_mean{k} = mean(yobs_all{k},2)' -4*pi*x{k}.^3/3;
        yobs_std{k}  = std(yobs_all{k},0,2)'; %)./sqrt(length(yobs_all{k}(1,:))
        yobs_ms{k}   = [yobs_mean{k}-yobs_std{k}; yobs_mean{k}+yobs_std{k}];
    end
end

ylo_min_all = cell2mat(ylo_min');
if(length(ylo_min_all(:,1))>1)
ylo_min_all = min(ylo_min_all);
end

yhi_max_all = cell2mat(yhi_max');
if(length(yhi_max_all(:,1))>1)
yhi_max_all = max(yhi_max_all);
end
ycsr_env_all  = [ylo_min_all;   yhi_max_all];


cmap=parula(6);
%figure('Name',name); 
hax=gca;

hold on
for k=1:length(A)

plot      (x{k},yobs_mean{k},'--','Color',cmap(k,:),'LineWidth',2)
plotshaded(x{k},yobs_ms{k},cmap(k,:)) %

end
plot      (x{1},ytheo{1},'--r','LineWidth',3)
plotshaded(x{1},ycsr_env_all ,[17 17 17]/255)

%plotshaded(x,yobs_env ,'b')
%xlabel('Neighborhood radius $r$ $(\mu m)$','Interpreter','latex')
%ylabel('$K(r)-\frac{4}{3}\pi r^3$','Interpreter','latex')
%legend({'Observed','Random','CSR envelope'},'Interpreter','latex','FontSize',14,...
%         'Location','northwest','NumColumns',1)
%hax.FontSize=16;
axis([0 1800 -1.6e+09 1.6e+09]) %1768

hold off
