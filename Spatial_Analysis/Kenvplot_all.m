%% Plots all the summary functions for three-dimensional point patterns.
% Author: Nikolaos M. Dimitriou, 
% McGill University, 2020

clear; clc; close all;

tp    = {'D0' 'D2' 'D5' 'D7' 'D9' 'D12' 'D14'};
% Series 1
group1 = {'*A*C*.csv','*A*E*.csv','*B*E*.csv','*B*N*.csv','*B*W*.csv','*F*W*.csv'};
gname1 = {'Control_s1_AC','Control_s1_AE','Control_s1_BE','Control_s1_BN','Control_s1_BW','Control_s1_FW'};

file1 = 'resEnv_sim_comp/K_env/';
gname = gname1;

%gname = [gname1s, gname1];

% Series 2
%

group2 = {%'*A*C*.csv','*A*N*.csv','*A*S*.csv','*A*W*.csv','*A*E*.csv',...
         %'*B*C*.csv','*B*N*.csv','*B*S*.csv','*B*W*.csv','*B*E*.csv',...
         %'*C*C*.csv','*C*N*.csv','*C*S*.csv','*C*W*.csv','*C*E*.csv',...
         %'*DD*C*.csv','*DD*N*.csv','*DD*S*.csv','*DD*W*.csv','*DD*E*.csv',...
                    '*E*N*.csv','*E*S*.csv','*E*W*.csv'           ,...
                    '*F*N*.csv',           '*F*W*.csv','*F*E*.csv'};
gname2 = {%'Pac_0p5_AC','Pac_0p5_AN','Pac_0p5_AS','Pac_0p5_AW','Pac_0p5_AE',...
          %'Pac_0p05_BC','Pac_0p05_BN','Pac_0p05_BS','Pac_0p05_BW','Pac_0p05_BE',...
          %'Pac_0p005_CC','Pac_0p005_CN','Pac_0p005_CS','Pac_0p005_CW','Pac_0p005_CE',...
          %'Pac_0p0005_DC','Pac_0p0005_DN','Pac_0p0005_DS','Pac_0p0005_DW','Pac_0p0005_DE',...
              'Control_s2_EN','Control_s2_ES','Control_s2_EW'     ,...
              'Control_s2_FN',     'Control_s2_FW','Control_s2_FE'};

          
gname = [gname2, gname1];

% File
file2 = 'resEnv_series_2_sim_comp/K_env/';
%


% CA series 1-2
file12s = 'resEnv_CA_s12/';

group1s = {'*ACs1*.csv','*AEs1*.csv','*BEs1*.csv','*BNs1*.csv','*BWs1*.csv','*FWs1*.csv'};
gname1s = {'CA_Control_s1_AC','CA_Control_s1_AE','CA_Control_s1_BE','CA_Control_s1_BN','CA_Control_s1_BW','CA_Control_s1_FW'};

%file1s = 'resEnv_CA_s1/run4_adhes-5_phenchange0/';


group2s = {'*ENs2*.csv','*ESs2*.csv','*EWs2*.csv','*FEs2*.csv','*FNs2*.csv','*FWs2*.csv'};
gname2s = {'CA_Control_s2_EN','CA_Control_s2_ES','CA_Control_s2_EW','CA_Control_s2_FE','CA_Control_s2_FN','CA_Control_s2_FW'};

%file2s = 'resEnv_CA_s2/K_env/';


gname_ca = [gname2s, gname1s];

run plotopt.m

figfile = '.'; %'Figures/CA_series_1/run4_adhes-5_phenchange0/';


%%
Kenv = struct;

Kenv = import_Kenv(file1,group1,tp,gname1,Kenv);
Kenv = import_Kenv(file2,group2,tp,gname2,Kenv);
Kenv = import_Kenv(file12s,group1s,tp,gname1s,Kenv);
Kenv = import_Kenv(file12s,group2s,tp,gname2s,Kenv);

%% Plot
c=1;
for i=1:length(group)
    h=figure('Name',gname{i},'Position', [100, 100, 1200, 900],'Visible','Off');
    ah = gobjects(7, 1);
    for j=1:length(tp)
       
        subplot(3,3,c);
        Kenvplot(Kenv.(tp{j}).(gname{i}));
        title(tp{j},'Interpreter','latex');
        c=c+1;
    end
    hh = mtit(' ');
    xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$');
    set(xlh, 'Visible', 'On');
    ylh=ylabel(hh.ah,'$K(r)-\frac{4}{3}\pi r^3$');
    set(ylh, 'Visible', 'On');
    hL=legend({'Observed','Random','CSR envelope'},'Interpreter','latex','FontSize',16,...
        'Location','southeastoutside','NumColumns',1);
    newPosition = [0.4 0.2 0.2 0.1];
    newUnits = 'normalized';
    set(hL,'Position', newPosition,'Units', newUnits);
    c=1;  
    
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(h,[figfile 'K_' gname{i} '.png'],'-dpng','-r350')
end
%% Envelope from all Paclitaxel 0.5uM distributions
h=figure('Name','Average K','Position', [100, 100, 1200, 900]);
c=1;
env_0p5 = {};
for i=1:length(tp)
   
    env{i}=struct2cell(Kenv.(tp{i}));
    % Paclitaxel 0.5um
    for j=1:5
        env_0p5{i}{j}=env{i}{j,1};
    end
    subplot(3,3,c);
    envelope(env_0p5{i});
    title(tp{i},'Interpreter','latex');
    c=c+1;
    
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On');
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On');
hL=legend({'$Pac\ 0.5 \mu M$','$SEM_{Pac\ 0.5 \mu M}$','Random','CSR envelope'},'Interpreter','latex','FontSize',16,...
        'Location','southeastoutside','NumColumns',1);
newPosition = [0.4 0.2 0.2 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,[figfile 'K_average_Pac_0p5.png'],'-dpng','-r350')

%% Envelope from all Paclitaxel 0.05uM distributions
h=figure('Name','Average K','Position', [100, 100, 1200, 900]);
c=1;
env_0p05 = {};
for i=1:length(tp)
    m=1;
    % Paclitaxel 0.05um
    for j=6:10
        env_0p05{i}{m}=env{i}{j,1};
        m=m+1;
    end
    subplot(3,3,c);
    envelope(env_0p05{i});
    title(tp{i},'Interpreter','latex');
    c=c+1;
    
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On');
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On');
hL=legend({'$Pac\ 0.05 \mu M$','$SEM_{Pac\ 0.05 \mu M}$','Random','CSR envelope'},'Interpreter','latex','FontSize',16,...
        'Location','southeastoutside','NumColumns',1);
newPosition = [0.4 0.2 0.2 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,[figfile 'K_average_Pac_0p05.png'],'-dpng','-r350')

%% Envelope from all Paclitaxel 0.005uM  distributions
h=figure('Name','Average K','Position', [100, 100, 1200, 900]);
c=1;
env_0p005 = {};
for i=1:length(tp)
    m=1;
    env{i}=struct2cell(Kenv.(tp{i}));
    % Paclitaxel 0.005um
    for j=11:15
        env_0p005{i}{m}=env{i}{j,1};
        m=m+1;
    end
    subplot(3,3,c);
    envelope(env_0p005{i});
    title(tp{i},'Interpreter','latex');
    c=c+1;
    
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On');
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On');
hL=legend({'$Pac\ 0.005 \mu M$','$SEM_{Pac\ 0.005 \mu M}$','Random','CSR envelope'},'Interpreter','latex','FontSize',16,...
        'Location','southeastoutside','NumColumns',1);
newPosition = [0.4 0.2 0.2 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,[figfile 'K_average_Pac_0p005.png'],'-dpng','-r350')

%% Envelope from all Paclitaxel 0.0005uM distributions
h=figure('Name','Average K','Position', [100, 100, 1200, 900]);
c=1;
env_0p0005 = {};
for i=1:length(tp)
    m=1;
    env{i}=struct2cell(Kenv.(tp{i}));
    % Paclitaxel 0.0005um
    for j=16:20
        env_0p0005{i}{m}=env{i}{j,1};
        m=m+1;
    end
    subplot(3,3,c);
    envelope(env_0p0005{i});
    title(tp{i},'Interpreter','latex');
    c=c+1;
    
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On');
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On');
hL=legend({'$Pac\ 0.0005 \mu M$','$SEM_{Pac\ 0.0005 \mu M}$','Random','CSR envelope'},'Interpreter','latex','FontSize',16,...
        'Location','southeastoutside','NumColumns',1);
newPosition = [0.4 0.2 0.2 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,[figfile 'K_average_Pac_0p0005.png'],'-dpng','-r350')

%% Envelope from all control distributions - series 1 and 2
run plotopt.m
h=figure('Name','Average K','Position', [100, 100, 1100, 900]);
c=1;
env_cont = {};
for i=1:length(tp)
    m=1;
    env{i}=struct2cell(Kenv.(tp{i}));
    % Control
    for j=1:12 %21:32
        env_cont{i}{m}=env{i}{j,1};
        m=m+1;
    end
    length(tp)
    subplot(3,3,c);
    envelope(env_cont{i});
    title(tp{i},'Interpreter','latex');
    c=c+1;
    
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(mm)$'); 
set(xlh, 'Visible', 'On');
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On');
hL=legend({'Observed','SEM$_{observed}$','Random','CSR envelope'},'Interpreter','latex','FontSize',24,...
        'Location','southeastoutside','NumColumns',1);
newPosition = [0.4 0.15 0.2 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h,['K_average_control_all.png'],'-dpng','-r350')


%% Envelope from control distributions CA series 1 and 2
run plotopt.m
%h=figure('Name','Average K','Position', [100, 100, 1000, 1000]);
c=1;
env_ca = {};
for i=1:length(tp)
    m=1;
    env{i}=struct2cell(Kenv.(tp{i}));
    % Control
    for j=13:24
        env_ca{i}{m}=env{i}{j,1};
        m=m+1;
    end

    %subplot(3,3,c);
    %envelope(env_cont{i});
    %title(tp{i},'Interpreter','latex');
    c=c+1;
    
end

%% Envelope from control distributions series 1 and 2
run plotopt.m
%h=figure('Name','Average K','Position', [100, 100, 1000, 1000]);
c=1;
env_s = {};
for i=1:length(tp)
    m=1;
    env{i}=struct2cell(Kenv.(tp{i}));
    % Control
    for j=1:12
        env_s{i}{m}=env{i}{j,1};
        m=m+1;
    end

    %subplot(3,3,c);
    %envelope(env_cont{i});
    %title(tp{i},'Interpreter','latex');
    c=c+1;
    
end

%% Plot K's average for each dataset
run plotopt.m
h=figure('Name','Average K','Position', [100, 100, 1000, 950]);
%h=figure('Name','Average K','Position', [100, 100, 800,700]);
c=1;

for i=1:length(tp)
    subplot(3,3,c);
    %envelopem({env_cont{i};env_0p0005{i};env_0p005{i};env_0p05{i};env_0p5{i}}); %;env_0p05{i};env_0p005{i};env_0p0005{i}
    envelopem({env_s{i};env_ca{i}});
    
    title(tp{i},'Interpreter','latex');
    c=c+1;
    
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On','FontSize',34);
ylh=ylabel(hh.ah,'$\langle K(r) \rangle - \frac{4}{3}\pi r^3$');
set(ylh, 'Visible', 'On','FontSize',34);
%hL=legend({'$Control$','$SEM_{Control}$','$Pac\ 0.0005 \mu M$','$SEM_{Pac\ 0.0005 \mu M}$','$Pac\ 0.005 \mu M$','$SEM_{Pac\ 0.005 \mu M}$','$Pac\ 0.05 \mu M$','$SEM_{Pac\ 0.05 \mu M}$','$Pac\ 0.5 \mu M$','$SEM_{Pac\ 0.5 \mu M}$',...
%    'Random','CSR envelope'},'Interpreter','latex','FontSize',20,...
%        'Location','southwest','NumColumns',2);
hL=legend({'Experiments','SD$_{Exp}$','Simulations','SD$_{Sim}$', 'Random','CSR envelope'},'Interpreter','latex','FontSize',20,...
            'Location','southeastoutside','NumColumns',1);     
            %'$SEM_{Control}$',...
            %'$Pac\ 0.5 \mu M$',...
            % '$SEM_{Pac\ 0.5 \mu M}$',...
            %'$Pac\ 0.05 \mu M$',... %'$SEM_{Pac\ 0.05 \mu M}$',...
            %'$Pac\ 0.005 \mu M$',... %'$SEM_{Pac\ 0.005 \mu M}$',...
            %'$Pac\ 0.0005 \mu M$','$SEM_{Pac\ 0.0005 \mu M}$',...
            
newPosition = [0.5 0.2 0.1 0.01];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(h,['K_average_all.png'],'-dpng','-r350')
print(h,['Figures/series12/K_average_all_sim_vs_exp_s12.png'],'-dpng','-r350')
%% Calculate differences of K between two subsequent days
%{
for i=1:length(tp)-1
    j=i+1;
    for k=1:length(group)
   
        dKenv.([tp{i} tp{j}]).(gname{k})=[Kenv.(tp{j}).(gname{k}).r, Kenv.(tp{j}).(gname{k}).obs - Kenv.(tp{i}).(gname{k}).obs];
        
    end
end

%% Plot the envelope of the differences
h=figure('Name','Average dK','Position', [100, 100, 1200, 900]);
for i=1:length(tp)-1
   
    j=i+1;
    denv{i}=struct2cell(dKenv.([tp{i} tp{j}]));
    subplot(3,3,i);
    diffenvelope(denv{i});
    title([tp{j} '-' tp{i}],'Interpreter','latex');
    
end
hh = mtit(' ');
xlh=xlabel(hh.ah,'Neighborhood radius $r$ $(\mu m)$'); 
set(xlh, 'Visible', 'On');
ylh=ylabel(hh.ah,'$\langle \Delta K(r) \rangle $');
set(ylh, 'Visible', 'On');
hL=legend({'Observed','Observed SEM'},'Interpreter','latex','FontSize',16,...
        'Location','southeastoutside','NumColumns',1);
%hL=legend(gname,'Interpreter','latex','FontSize',16,...
%        'Location','southeastoutside','NumColumns',1);
newPosition = [0.8 0.2 0.1 0.1];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h,['Figures/diffK_average.png'],'-dpng','-r0')
%}
