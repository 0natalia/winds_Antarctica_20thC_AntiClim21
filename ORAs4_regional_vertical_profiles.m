%__________________________________________________________________________
% 		    THETA TEMPERATURE ANALYSIS IN THE ANTARCTIC COAST             %%
%               Oceanic reanalysis data from ECMWF/ORAs4                   %
% Monthly, 42 depth levels, VAR = THETA0, DEPTH, MASK, TIME, LON, LAT %

% Nat?lia Silva; natalia3.silva@usp.br
% 2020
%__________________________________________________________________________
 
clear all; close all

dado = ('thetao_oras4_1m_1958_grid_1x1.nc');b % choose any to create mask
lat = ncread(dado,'lat'); l = find(lat<-59); lat = lat(l);
lon = ncread(dado,'lon');
time = ncread(dado,'time'); tempo = datenum(datetime(1957,9,1)+days(time)); clear time
dep = ncread(dado,'depth'); dp = find(dep < 1000); dep = dep(dp); % 1000m
mask = ncread(dado,'mask'); % land = 0; oc = 1;
clear dado

% Read data
d = dir('db1/reanalise/ocn/ORAS4/THETA/');
t = [];
for i = 3:length(d)
    ora = d(i).name;
    theta = ncread(ora,'thetao'); theta = theta(:,l,dp,:);
    theta = squeeze(mean(theta,4));
    t = cat(4,t,theta);
end
clear i; clear ora; clear l; clear d

% mean Before75 and After75
t_B75 = mean(t(:,:,:,1:18),4); 
t_A75 = mean(t(:,:,:,18:end),4);
v=(-2:0.1:2); vv = (-0.5:0.01:0.5);

%% %%%%%%% 120E - East
ll = find(lon ==120.5); la = find(lat>-68);

% for i = 1:6
t_reg_B75 = squeeze(t_B75(ll,la,:));
t_reg_A75 = squeeze(t_A75(ll,la,:));
dif = t_reg_A75 - t_reg_B75;

figure('color',[1 1 1],'position',[108 305 1300 600]);
%     subplot(1,3,1)
%     contourf(lat(la),-dep,t_reg_B75',v); hold on
%     contour(lat(la),-dep,t_reg_B75',v)
%     colormap(flipud(cbrewer('div','Spectral',40)))
%     title('\Theta_{B75} 300 < lon < 360'); ylabel('z (m)','Fontsize',9); xlabel('lat')
%     caxis([-2 2])
%     colorbar;title(colorbar,'^oC')
% 
%     subplot(1,3,2);
%     contourf(lat(la),-dep,t_reg_A75',v); hold on
%     contour(lat(la),-dep,t_reg_A75',v)
%     title('\Theta_{A75}  300 < lon < 360'); ylabel('z (m)','Fontsize',9); xlabel('lat')
%     caxis([-2 2])
%     colorbar;title(colorbar,'^oC')
    
subplot(1,3,3)
contourf(lat(la),-dep,dif',vv); hold on
contour(lat(la),-dep,dif',vv)
colormap(cmocean('balance'))
title('\Delta\Theta 120^oE'); ylabel('z (m)','Fontsize',9); xlabel('lat')
caxis([-0.5 0.5])
colorbar;title(colorbar,'^oC')

%     ll = ll+60;
% end
%clear i

%%%%%%%%%%%%%%%% 112 W

ll = find(lon ==248.5); la = find(lat>-75);

t_reg_B75 = squeeze(t_B75(ll,la,:));
t_reg_A75 = squeeze(t_A75(ll,la,:));
dif = t_reg_A75 - t_reg_B75;

figure('color',[1 1 1],'position',[108 305 1300 600])
subplot(1,3,3)
contourf(lat(la),-dep,dif',vv); hold on
contour(lat(la),-dep,dif',vv)
colormap(cmocean('balance'))
title('\Delta\Theta 112^oW'); ylabel('z (m)','Fontsize',9); xlabel('lat')
caxis([-0.5 0.5])
colorbar;title(colorbar,'^oC')

%%%%%%%%%%%%%%%% 68 S

ll = find(lon>270 & lon<310); la = find(lat==-68.5);

t_reg_B75 = squeeze(t_B75(ll,la,:));
t_reg_A75 = squeeze(t_A75(ll,la,:));
dif = t_reg_A75 - t_reg_B75;

figure('color',[1 1 1],'position',[108 305 1300 600])
subplot(1,3,3)
contourf(lon(ll),-dep,dif',vv); hold on
contour(lon(ll),-dep,dif',vv)
colormap(cmocean('balance'))
title('\Delta\Theta 68^oS'); ylabel('z (m)','Fontsize',9); xlabel('lon')
caxis([-0.5 0.5])
colorbar;title(colorbar,'^oC')

