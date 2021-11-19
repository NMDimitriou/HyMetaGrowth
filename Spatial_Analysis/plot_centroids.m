%% 3D plot of centroids
clear; clc;

tp     = {'D0' 'D2' 'D5' 'D7' 'D9' 'D12' 'D14'};
ts     = {'D0' 'D1' 'D2' 'D3' 'D4' 'D5' 'D6' 'D7' 'D8' 'D9' 'D10' 'D11' 'D12' 'D13' 'D14'};
% Series 1
%
%group = {'A*-C.tif_plot3D.fig','A*-E.tif_plot3D.fig','B*-E.tif_plot3D.fig',...
%    'B*-N.tif_plot3D.fig','B*-W.tif_plot3D.fig','F*-W.tif_plot3D.fig'};

% Series 1
group1 = {'AD*C.tif.txt','AD*E.tif.txt','BD*E.tif.txt','BD*N.tif.txt','BD*W.tif.txt','FD*W.tif.txt'};

groups1= {'num_cells_AC.txt','num_cells_AE.txt','num_cells_BE.txt','num_cells_BN.txt','num_cells_BW.txt','num_cells_FW.txt'};

groupc1= {'ca_res_AC.csv','ca_res_AE.csv','ca_res_BE.csv','ca_res_BN.csv','ca_res_BW.csv','ca_res_FW.csv'};


gname1  = {'Control_ACs1','Control_AEs1','Control_BEs1','Control_BNs1','Control_BWs1','Control_FWs1'};

       
%gname1  = {'Sample 1','Sample 2','Sample 3','Sample 4','Sample 5','Sample 6',...
%           'Sample 7','Sample 8','Sample 9','Sample 10','Sample 11','Sample 12'};

% CA series 1,2
%
group_sim = 'run4_adhes-5_phenchange0_set12_best_params/';


% Series 2

group2 = {%'A*C*.txt','A*N*.txt','A*S*.txt','A*W*.txt','A*E*.txt',...
         %'B*C*.txt','B*N*.txt','B*S*.txt','B*W*.txt','B*E*.txt',...
         %'C*C*.txt','C*N*.txt','C*S*.txt','C*W*.txt','C*E*.txt',...
         %'D*C*.txt','D*N*.txt','D*S*.txt','D*W*.txt','D*E*.txt',...
                    'E*N*.txt','E*S*.txt','E*W*.txt'           ,...
                    'F*N*.txt',           'F*W*.txt','F*E*.txt'};

groups2 = {'num_cells_Control_s2_EN.txt','num_cells_Control_s2_ES.txt','num_cells_Control_s2_EW.txt',...
    'num_cells_Control_s2_FN.txt','num_cells_Control_s2_FW.txt','num_cells_Control_s2_FE.txt'};                
                
groupc2 = {'ca_res_Control_s2_EN.csv','ca_res_Control_s2_ES.csv','ca_res_Control_s2_EW.csv',...
           'ca_res_Control_s2_FN.csv','ca_res_Control_s2_FW.csv','ca_res_Control_s2_FE.csv'};                
                
gname2 = {%'Pac_0p5_AC','Pac_0p5_AN','Pac_0p5_AS','Pac_0p5_AW','Pac_0p5_AE',...
          %'Pac_0p05_BC','Pac_0p05_BN','Pac_0p05_BS','Pac_0p05_BW','Pac_0p05_BE',...
         %'Pac_0p005_CC','Pac_0p005_CN','Pac_0p005_CS','Pac_0p005_CW','Pac_0p005_CE',...
         %'Pac_0p0005_DC','Pac_0p0005_DN','Pac_0p0005_DS','Pac_0p0005_DW','Pac_0p0005_DE',...
              'Control_ENs2','Control_ESs2','Control_EWs2'     ,...
              'Control_FNs2',     'Control_FWs2','Control_FEs2'};

time  = [0 2 5 7 9 12 14];

run plotopt.m

file1e  = 'res_coord_scaled/';
file2e  = 'res_coord_series_2_scaled/';
file_s  = 'res_coord_sim_series_12/run4_adhes-5_phenchange0/';

ntp     = length(tp);

%%

for i=1:length(group1)
   
    for j=1:length(tp)
       
        samp1{i,j}=dir(['res_coord_scaled/' tp{j} '/' group1{i}]);
        coord1.(gname1{i}).(tp{j})=readmatrix([samp1{i,j}.folder '/' samp1{i,j}.name]);
        count1.(gname1{i}).(tp{j})=length(coord1.(gname1{i}).(tp{j})(:,1));

    end
    samps1{i}=dir([group_sim groups1{i}]);
    sampc1{i}=dir([group_sim groupc1{i}]);
    ca_count1.(gname1{i}) = readmatrix([samps1{i}.folder '/'  samps1{i}.name]);
    ca_coord1.(gname1{i}) = readmatrix([sampc1{i}.folder '/'  sampc1{i}.name]);
    
end

for i=1:length(group2)  
    for j=1:length(tp)
       
        samp2{i,j}=dir(['res_coord_series_2_scaled/' tp{j} '/' group2{i}]);
        coord2.(gname2{i}).(tp{j})=readmatrix([samp2{i,j}.folder '/' samp2{i,j}.name]);
        count2.(gname2{i}).(tp{j})=length(coord2.(gname2{i}).(tp{j})(:,1));

    end
    if(i<7)
        samps2{i}=dir([group_sim groups2{i}]);
        sampc2{i}=dir([group_sim groupc2{i}]);
        ca_count2.(gname2{i}) = readmatrix([samps2{i}.folder '/'  samps2{i}.name]);
        ca_coord2.(gname2{i}) = readmatrix([sampc2{i}.folder '/'  sampc2{i}.name]);
    end
end

%%
% Split time-points of CA data
ca_coord1_tsp = split_tp(groupc1, gname1, ca_coord1, ts);
ca_coord2_tsp = split_tp(groupc2, gname2, ca_coord2, ts);

%%

% Reformat coordinates
sz=[480,480,176];
h=2500/480; %scale factor

for i=1:length(ts)
    mkdir([group_sim 'scale_coord/' ts{i} '/']);
end

%

ca_coord1_tsp = save_sim_split_coord(group_sim,groupc1,ca_coord1_tsp,gname1,ts,h,sz);
ca_coord2_tsp = save_sim_split_coord(group_sim,groupc2,ca_coord2_tsp,gname2,ts,h,sz);


%for i=1:length(group1)
%    for j=1:length(tp)
        
%        A(j)=hgload([samp1{i,j}.folder '/' samp1{i,j}.name]);
        
        
%    end
%end
%% Merge controls and count all
gc1=struct2array(count1);
gc2=struct2array(count2);
gc1=struct2table(gc1);
gc2=struct2table(gc2);
gc1=table2array(gc1);
gc2=table2array(gc2);

gc_control = gc1;
%gc_control = [gc1; gc2];
%gc_control = [gc1; gc2(21:26, :)];
%gc_pac_0p5 = gc2(1:5,:);
%gc_pac_0p05 = gc2(6:10,:);
%gc_pac_0p005 = gc2(11:15,:);
%gc_pac_0p0005 = gc2(16:20,:);

gcm_control   = mean(gc_control);
gcstd_control = std(gc_control);

%gcm_pac_0p5 = mean(gc_pac_0p5);
%gcm_pac_0p05 = mean(gc_pac_0p05);
%gcm_pac_0p005 = mean(gc_pac_0p005);
%gcm_pac_0p0005 = mean(gc_pac_0p0005);

%gcstd_pac_0p5 = std(gc_pac_0p5);
%gcstd_pac_0p05 = std(gc_pac_0p05);
%gcstd_pac_0p005 = std(gc_pac_0p005);
%gcstd_pac_0p0005 = std(gc_pac_0p0005);

% find std of ca simulated data

ts = (0:0.001:14);

% interpolate accoring to ts
%
ca_count_table = [];
ca_count_table = merge_counts(ca_count_table,groups1,ca_count1,gname1,ts,0);

%
%ca_count_table = merge_counts(ca_count_table,groups2,ca_count2,gname2,ts,6);
%

ca_count_mean= mean(ca_count_table,2);
ca_count_std = std(ca_count_table,0,2);

%% Plot growth
f1=figure('Name','Growth','Position', [100, 100, 800, 700]);
hold on
h1=plot(ts,ca_count_mean,'b','LineWidth',4,'DisplayName','Simulations');
plotshaded(ts,[ca_count_mean'-ca_count_std';ca_count_mean'+ca_count_std'],'b');
h2=errorbar(time,gcm_control,gcstd_control,'.-r','LineWidth',8,'MarkerSize',35,'DisplayName','Experiments');
legend([h1 h2],'Location','northwest','FontSize',36)
%errorbar(time,gcm_pac_0p0005,gcstd_pac_0p0005,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.0005\ \mu$M')
%errorbar(time,gcm_pac_0p005,gcstd_pac_0p005,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.005\ \mu$M')
%errorbar(time,gcm_pac_0p05,gcstd_pac_0p05,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.05\ \mu$M')
%errorbar(time,gcm_pac_0p5,gcstd_pac_0p5,'.-','LineWidth',6,'MarkerSize',30,'DisplayName','$0.5\ \mu$M')
xlabel('time $(days)$','Interpreter','latex')
ylabel('Nuclei count','Interpreter','latex')
%legend('Location','NorthWest')
axis([0 15 -inf inf])
hold off
set(f1,'Units','Inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
disp('Saving...')
print(f1,['nuclei_count_s1.png'],'-r350','-dpng')

%% Plot centroids Series 1-2 experiment
plot_centr(group1,gname1,coord1,tp,'centroids');
plot_centr(group2,gname2,coord2,tp,'centroids');


clf; close all

%% Plot centroids Series 1-2 CA
plot_centr_CA(groupc1,gname1,ca_coord1_tsp,tp,'centroids_CAs1');
plot_centr_CA(groupc2,{gname2{21:end}},ca_coord2_tsp,tp,'centroids_CAs2');


clf; close all

%% Plot CA and experiment together
run plotopt.m
%plot_centr_CA_exp(groupc1,ca_coord1_tsp,coord1,gname1,tp,'centroids_CA_exp_s1')
plot_centr_CA_exp(groupc2,ca_coord2_tsp,coord2,gname2,tp,'centroids_CA_exp_s2')


clf; close all