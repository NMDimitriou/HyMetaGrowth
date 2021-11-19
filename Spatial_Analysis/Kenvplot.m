function varargout = Kenvplot(A)
%% Plots the summary function for a three-dimensional point pattern.
% input: the matrix obtained from envelope.pp3 function of spatstat package
% of R
% output: the plot
% Author: Nikolaos M. Dimitriou, 
% McGill University, 2020

x=A.r';
yobs=A.obs';
yobs=yobs-4*pi*x.^3/3;
ytheo=A.theo';
ytheo=ytheo-4*pi*x.^3/3;
yenv=[A.lo, A.hi];
yenv=yenv';
yenv=yenv-4*pi*x.^3/3;

%figure('Name',name); 
hold on
hax=gca;
plot      (x,yobs ,'LineWidth',4)
plot      (x,ytheo,'--','LineWidth',3)
plotshaded(x,yenv ,[17 17 17]/255)
%xlabel('Neighborhood radius $r$ $(\mu m)$','Interpreter','latex')
%ylabel('$K(r)-\frac{4}{3}\pi r^3$','Interpreter','latex')
%legend({'Observed','Random','CSR envelope'},'Interpreter','latex','FontSize',14,...
%         'Location','northwest','NumColumns',1)

hax.FontSize=16;
axis([-inf inf -2.5e+09 2.5e+09])
hold off
