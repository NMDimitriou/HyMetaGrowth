%% Calculate and plot distances of nuclei for all time points, in each sample
% Author: Nikolaos M. Dimitriou, 
% McGill University, 2020
clear; clc; close all;

% the parent directory that contains the subdirectories of the time-points
file1e  = 'res_coord_scaled/';
%file2e  = 'res_coord_series_2_scaled/';
file_s  = 'res_coord_sim_series_12/run4_adhes-5_phenchange0_set12_best_params/scale_coord/';

% the list of the time-point subdirectories that contain the data
tp     = {'D0' 'D2' 'D5' 'D7' 'D9' 'D12' 'D14'};
time   = [0 2 5 7 9 12 14];

% the name of the data files for the coordinates
% Series 1
group1 = {'A*C*.txt','A*E*.txt','B*E*.txt','B*N*.txt','B*W*.txt','F*W*.txt'};

% Series 2
group2 = {%'A*C*.txt','A*N*.txt','A*S*.txt','A*W*.txt','A*E*.txt',...
          %'B*C*.txt','B*N*.txt','B*S*.txt','B*W*.txt','B*E*.txt',...
          %'C*C*.txt','C*N*.txt','C*S*.txt','C*W*.txt','C*E*.txt',...
          %'D*C*.txt','D*N*.txt','D*S*.txt','D*W*.txt','D*E*.txt',...
                    'E*N*.txt','E*S*.txt','E*W*.txt'           ,...
                    'F*N*.txt',           'F*W*.txt','F*E*.txt'};


groupcs1={'CA_coord_Control_ACs1*.txt','CA_coord_Control_AEs1*.txt','CA_coord_Control_BEs1*.txt',...
          'CA_coord_Control_BNs1*.txt','CA_coord_Control_BWs1*.txt','CA_coord_Control_FWs1*.txt'};
        
groupcs2={'CA_coord_Control_ENs2*.txt','CA_coord_Control_ESs2*.txt','CA_coord_Control_EWs2*.txt',...
          'CA_coord_Control_FNs2*.txt','CA_coord_Control_FWs2*.txt','CA_coord_Control_FEs2*.txt'};
      
gname1 = {'Sample_1', 'Sample_2', 'Sample_3', 'Sample_4' , 'Sample_5' , 'Sample_6' };
gname2 = {'Sample_7', 'Sample_8', 'Sample_9', 'Sample_10', 'Sample_11', 'Sample_12' };

gname_ca1 = {'CA_Sample_1', 'CA_Sample_2', 'CA_Sample_3', 'CA_Sample_4', 'CA_Sample_5', 'CA_Sample_6' };
gname_ca2 = {'CA_Sample_7', 'CA_Sample_8', 'CA_Sample_9', 'CA_Sample_10', 'CA_Sample_11', 'CA_Sample_12' };

gname = [gname1, gname2, gname_ca1, gname_ca2];
datgroup={'Experiments','Simulations'};

nsamp1  = [length(group1), length(groupcs1)];
nsamp2  = [length(group2), length(groupcs2)];
ntp    = length(tp);

vname = {'ACs1','AEs1','BEs1','BNs1','BWs1','FWs1',...
         'ENs2','ESs2','EWs2','FNs2','FWs2','FEs2'};

run plotopt.m

%% Import coordinates
disp('1. Importing coordinates...')
% Series 1-2
samp1    = cell(length(vname),7);
sampc1   = cell(length(vname),7);
coord    = struct;
count    = struct;
ca_coord = struct;
[samp1,sampc1,coord,count,ca_coord] = import_coord(nsamp1(1),tp,ntp,file1e,file_s,group1,gname1,groupcs1,gname_ca1,...
    samp1,sampc1,coord,count,ca_coord,0);
[samp1,sampc1,coord,count,ca_coord] = import_coord(nsamp2(2),tp,ntp,file2e,file_s,group2,gname2,groupcs2,gname_ca2,...
    samp1,sampc1,coord,count,ca_coord,6);

disp('Finished importing coordinates...')

%% Plot histogram of cells across z
disp('2. Plotting histograms...')
run plotopt.m

plot_hist_z(coord,ca_coord,tp,'hist_z_s12_exp_ca')

disp('4. Finished plotting histograms...')

%% Calculate distances and save them
disp('4. Calculating and saving distances...')
gnam    = [gname1;gname2];
gnam_ca = [gname_ca1;gname_ca2];
%parpool(nsamp(1)); %memory is a thing here
for i=1:(nsamp1(1)+nsamp2(2))
    for k=1:ntp
   
        calc_ind_knd(   coord.(gnam{i}).(tp{k}), samp1{i,k}.name)
        calc_ind_knd(ca_coord.(gnam_ca{i}).(tp{k}),sampc1{i,k}.name)
    
    end
end
disp('Finished calculating and saving distances!')


%% Estimate variations between distributions
% import distances
disp('5. Importing distances...')
for i=1:nsamp1(1) + nsamp2(2)
    for k=1:ntp

    sampINDist{i,k}=dir(['Distances/INDist/' 'INDist_' samp1{i,k}.name '.bin' ]);
    sampKNDist{i,k}=dir(['Distances/KNDist/' 'KNDist_' samp1{i,k}.name '.bin']);

    fileID=fopen([sampINDist{i,k}.folder '/' sampINDist{i,k}.name],'r');
    INDist.(gnam{i}).(tp{k})=fread(fileID,'double');
    fclose(fileID);

    fileID=fopen([sampKNDist{i,k}.folder '/' sampKNDist{i,k}.name],'r');
    KNDist.(gnam{i}).(tp{k})=fread(fileID,'double');
    fclose(fileID);

    end
end

gnam_ca = [gname_ca1;gname_ca2];
for i=1:nsamp1(1) + nsamp2(2)
    for k=1:ntp

    sampINDistCA{i,k}=dir(['Distances/INDist/' 'INDist_' sampc1{i,k}.name '.bin' ]);
    sampKNDistCA{i,k}=dir(['Distances/KNDist/' 'KNDist_' sampc1{i,k}.name '.bin']);

    fileID=fopen([sampINDistCA{i,k}.folder '/' sampINDistCA{i,k}.name],'r');
    INDist.(gnam_ca{i}).(tp{k})=fread(fileID,'double');
    fclose(fileID);

    fileID=fopen([sampKNDistCA{i,k}.folder '/' sampKNDistCA{i,k}.name],'r');
    KNDist.(gnam_ca{i}).(tp{k})=fread(fileID,'double');
    fclose(fileID);

    end
end

disp('Finished importing distances!')

%% calculate the cosine similarity between distributions of simulations and experiments

disp('6. Estimating similarity sim vs exp...')
calc_dist_var_sim_exp(INDist,KNDist,nsamp1(1)+nsamp2(2),tp,vname);
disp('Finished estimating similarity sim vs exp!')


