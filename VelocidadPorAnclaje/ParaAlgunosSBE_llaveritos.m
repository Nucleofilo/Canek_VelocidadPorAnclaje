clr; 
% algunos sensores de temperautra no tienen presion y estan en anclajes
% que realmente no se mueven pero si registran cambios de presion
%
load T:\Chamba\ReconstruccionTemperatura\Redundancy\Offset\Gradients\Clim_emp_dTdp.mat

%

load I:\TesisDoc\datos\CNK_extractions\P3\A\CNK57_YUC3 
{anchor.name}';

tr = smoothdata(anchor(1).T, 'movmean', 2*24);
tim = anchor(1).time;
tc = smoothdata(anchor(2).T, 'movmean', 2*24);
pr = smoothdata(anchor(1).p, 'movmean', 2*24);
pc = smoothdata(anchor(2).p, 'movmean', 2*24);

xr = anchor(1).lola(1);

dt = nanmean(tr(500:6000)) - nanmean(tc(500:6000));
Tm = nanmean(tr(:));

dTdpr = interp2(elo, eT, dTdp_c, xr,Tm );

dz = dt/abs(dTdpr);


anchor(2).p = anchor(1).p - dz

%
save('I:\TesisDoc\datos\CNK_extractions\P3\A\CNK57_YUC3.mat', 'anchor');



%%


clr; 
% algunos sensores de temperautra no tienen presion y estan en anclajes
% que realmente no se mueven pero si registran cambios de presion
%
load T:\Chamba\ReconstruccionTemperatura\Redundancy\Offset\Gradients\Clim_emp_dTdp.mat
%
load I:\TesisDoc\datos\CNK_extractions\P3\A\CNK42_YUC1 
{anchor.name}';

tr = smoothdata(anchor(1).T, 'movmean', 2*24);
tim = anchor(1).time;
tc = smoothdata(anchor(2).T, 'movmean', 2*24);
pr = smoothdata(anchor(1).p, 'movmean', 2*24);
pc = smoothdata(anchor(2).p, 'movmean', 2*24);
xr = anchor(1).lola(1);
dt = nanmean(tr(500:6000)) - nanmean(tc(500:6000));
Tm = nanmean(tr(:));
dTdpr = interp2(elo, eT, dTdp_c, xr+0.35,Tm );
dz = dt/abs(dTdpr);
anchor(2).p = anchor(1).p - dz;
%%
save('I:\TesisDoc\datos\CNK_extractions\P3\A\CNK42_YUC1.mat', 'anchor');



%%
close all;
plot(anchor(1).p - dz)
hold on;
plot(anchor(1).p)


%