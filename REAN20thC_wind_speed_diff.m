%__________________________________________________________________________
%                       AntClim21 Workshop - BAS                          %
%               Wind changes in Antarctica during the 20th C              %
%              Reanalysis data from NOAA-20CR e ECMWF-ERA20C              %
%                            READER obs data                              %

% This script reads wind and SLP data; calculates the wind speed mean for
% the time series' 20 first and last years, and its difference; Wind speed
% differences for the entire austral region, and specfic regions (Wilkes
% Glaciers, Amundsen Sea and Western Antarctic Peninsula)

% Natalia Silva; natalia3.silva@usp.br
% 2020
%__________________________________________________________________________

clear all; close all

addpath('/home/natalia/rotinas_MATLAB')
addpath('/home/natalia/rotinas_MATLAB/m_map')
addpath('home/natalia/READER_Data/Wind_Speed')
addpath('home/natalia/rotinas_MATLAB/CBrewer_colormap/')

%% ECMWF- ERA20C
clear all; close all

era20u = ('home/natalia/reanalises/ERA20C/era_uwnd.nc');
era20v = ('home/natalia/reanalises/ERA20C/era_vwnd.nc');
pera = ('home/natalia/reanalises/ERA20C/era_slp.nc');

lat = ncread(era20u,'lat'); l = find(lat<-54); lat = lat(l);
lon = ncread(era20v,'lon'); 
time2 = ncread(era20u,'time'); tempo = datenum(datetime(1900,1,1)+hours(time2)); clear time2 

u2 = ncread(era20u,'u'); v2 = ncread(era20v,'v'); % nivel 5 = superficie
u2 = squeeze(u2(:,l,5,:)); v2 = squeeze(v2(:,l,5,:)); wsera = hypot(u2,v2);

u2 = [u2(:,1,:) u2]; v2 = [v2(:,1,:) v2];
wsera = [wsera(:,1,:) wsera]; lat = [-90; lat];
u2 = [u2; u2(end,:,:)]; v2 = [v2; v2(end,:,:)];
wsera = [wsera; wsera(end,:,:)]; lon = [lon; 180];

[wsera_anual,~] = downsample_ts(wsera,tempo,'year');
[u2_anual,~] = downsample_ts(u2,tempo,'year');
[v2_anual,~] = downsample_ts(v2,tempo,'year');
clear era20u; clear era20v 

% NOAA - CIRES 20CR
noaav = ('home/natalia/reanalises/NOAA20CR/vwnd.10m.mon.mean.nc');
noaau = ('home/natalia/reanalises/NOAA20CR/uwnd.10m.mon.mean.nc');
pnoaa = ('home/natalia/reanalises/NOAA20CR/slp.mon.mean.nc');

v1 = ncread(noaav,'vwnd'); v1 = v1(:,l,349:1680);
u1 = ncread(noaau,'uwnd'); u1 = u1(:,l,349:1680);  clear l;
wsnoaa = hypot(u1,v1);

u1 = [u1(:,1,:) u1]; v1 = [v1(:,1,:) v1];
wsnoaa = [wsnoaa(:,1,:) wsnoaa];
u1 = [u1; u1(end,:,:)]; v1 = [v1; v1(end,:,:)];
wsnoaa = [wsnoaa; wsnoaa(end,:,:)];

[wsnoaa_anual,tanos] = downsample_ts(wsnoaa,tempo,'year');
[u1_anual,~] = downsample_ts(u1,tempo,'year');
[v1_anual,~] = downsample_ts(v1,tempo,'year');
clear noaav; clear noaau; %clear tempo

% READER
d = dir('/home/natalia/READER_Data/Wind_Speed/W*.txt');

ws_reader = [];
for i = 1:length(d)    
	est = load(d(i).name); anual = [];  %disp(d(i).name)
    fim = find(est(:,1)==2010); 
    for j = 1:length(1:fim);    
        anual = [anual nanmean(est(j,2:13))];      
    end
    clear j
    ws_reader(end+1) = nanmean(anual);
end
clear i; clear anual; clear d; clear est; clear fim

ws_reader = ws_reader*0.5144; % converter de knots -> m/s

dado_coord = ('/home/natalia/READER_Data/coordenadas_estacoesREADER_vento.xlsx');
[~,latr,~] = xlsread(dado_coord,1,'B2:B20'); latr = str2double(latr);
[~,lonr,~] = xlsread(dado_coord,1,'C2:C20'); lonr = str2double(lonr);
clear dado_coord

clear pera; clear pnoaa
%% NOAA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wsmed1 = mean(wsnoaa(:,:,601:end),3);
u1med = mean(u1,3); v1med = mean(v1,3);

wsnoaa_p20 = mean(wsnoaa_anual(:,:,1:20),3); wsnoaa_u20 = mean(wsnoaa_anual(:,:,end-19:end),3);
dif_wsnoaa = wsnoaa_u20 - wsnoaa_p20;
u1_p20 = mean(u1_anual(:,:,1:20),3); u1_u20 = mean(u1_anual(:,:,end-19:end),3);
dif_u1 = u1_u20 - u1_p20;
v1_p20 = mean(v1_anual(:,:,1:20),3); v1_u20 = mean(v1_anual(:,:,end-19:end),3);
dif_v1 = v1_u20 - v1_p20;
clear wsnoaa_p20; clear wsnoaa_u20; clear u1_anual; clear v1_anual
clear u1_p20; clear u1_u20; clear v1_p20; clear v1_u20;
clear wsnoaa_anual

% ERA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wsmed = mean(wsera(:,:,601:end),3);
u2med = mean(u2,3); v2med = mean(v2,3);

wsera_p20 = mean(wsera_anual(:,:,1:20),3); wsera_u20 = mean(wsera_anual(:,:,end-19:end),3);
dif_wsera = wsera_u20 - wsera_p20;
u2_p20 = mean(u2_anual(:,:,1:20),3); u2_u20 = mean(u2_anual(:,:,end-19:end),3);
dif_u2 = u2_u20 - u2_p20;
v2_p20 = mean(v2_anual(:,:,1:20),3); v2_u20 = mean(v2_anual(:,:,end-19:end),3);
dif_v2 = v2_u20 - v2_p20;
clear wsera_p20; clear wsera_u20; clear u2_p20; clear u2_u20; 
clear v2_p20; clear v2_u20; clear u2_anual; clear v2_anual;
clear wsera_anual
 
[x,y] = meshgrid(lon,lat);
ttt=(0.8:0.1:12); pp=(-2.5:0.05:2.5);

%% MAPAS

figure('color',[1 1 1],'position',[10 305 940 815])
m_proj('stereographic','lat',-90,'long',0,'radius',35); 
colormap(cmocean('dense'))
m_contourf(lon,lat,wsmed1',ttt); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s')   
m_contour(lon,lat,wsmed1',ttt)  
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),u1med(1:2:end,1:2:end)',...
    v1med(1:2:end,1:2:end)',0.4,'color',[0.4 0.4 0.4]);
m_coast('color','k','linewidth',1.5); caxis([0.8 12]); 
m_scatter(lonr,latr,ws_reader*30, ws_reader, 'filled', 'MarkerEdgeColor', [1 1 1], 'LineWidth',0.9)
m_grid('xtick',[0 45 90 135 180 -45 -90 -135],'tickdir','out','ytick',...
    [-60 -75 -90],'tickdir','out', 'xaxisLocation', 'top', ...
    'yaxisLocation', 'middle','fontsize',7,'linestyle',':','linewidth',0.5);
legend(p, '12.5 m/s', 'location', 'southwest')
title('Wind Speed (1950:2010) - NOAA 20CR / READER','Position',[0 0.7 1])
clear t; clear tt;

figure('color',[1 1 1],'position',[10 305 940 815]);
m_proj('stereographic','lat',-90,'long',0,'radius',35);   
colormap(m_colmap('diverging'))
m_contourf(lon,lat,dif_wsnoaa',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsnoaa',pp)  
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u1(1:2:end,1:2:end)',...
    dif_v1(1:2:end,1:2:end)',0.4,'color',[0.5 0.5 0.5]);
m_grid('xtick',[0 45 90 135 180 -45 -90 -135],'tickdir','out','ytick',...
    [-60 -75 -90],'tickdir','out', 'xaxisLocation', 'top', ...
    'yaxisLocation', 'middle','fontsize',7,'linestyle',':','linewidth',0.5);
legend(p, '12.5 m/s', 'Location', 'southwest')
title('\DeltaWS (20L-20F) - NOAA 20CR','Position',[0 0.7 1]);
caxis([-2.5 2.5]); 
clear t; clear tt; clear ans;
clear u1med; clear v1med; clear wsmed1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('color',[1 1 1],'position',[10 305 940 815]);
m_proj('stereographic','lat',-90,'long',0,'radius',35);
colormap(cmocean('dense'))
m_contourf(lon,lat,wsmed',ttt); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s')   
m_contour(lon,lat,wsmed',ttt); hold on   
m_coast('color','k','linewidth',1.5); caxis([0.8 12]); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),u2med(1:2:end,1:2:end)',...
    v2med(1:2:end,1:2:end)',0.4,'color',[0.4 0.4 0.4]);
m_scatter(lonr,latr,ws_reader*30, ws_reader, 'filled', 'MarkerEdgeColor', [1 1 1], 'LineWidth',0.9)
m_grid('xtick',[0 45 90 135 180 -45 -90 -135],'tickdir','out','ytick',...
    [-60 -75 -90],'tickdir','out', 'xaxisLocation', 'top', ...
    'yaxisLocation', 'middle','fontsize',7,'linestyle',':','linewidth',0.5);
legend(p, '12.5 m/s', 'Location', 'southwest')
title('Wind Speed (1950:2010) - ERA 20C / READER','Position',[0 0.7 1])

clear t; clear tt; clear c; %clear ws_reader

figure('color',[1 1 1],'position',[10 305 940 815]);
m_proj('stereographic','lat',-90,'long',0,'radius',35);   
colormap(m_colmap('diverging'))
m_contourf(lon,lat,dif_wsera',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsera',pp)
m_coast('color','k','linewidth',1.5); caxis([-2.5 2.5]); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u2(1:2:end,1:2:end)',...
    dif_v2(1:2:end,1:2:end)',0.4,'color',[0.5 0.5 0.5]);
m_grid('xtick',[0 45 90 135 180 -45 -90 -135],'tickdir','out','ytick',...
    [-60 -75 -90],'tickdir','out', 'xaxisLocation', 'top', ...
    'yaxisLocation', 'middle','fontsize',7,'linestyle',':','linewidth',0.5);
legend(p, '12.5 m/s', 'Location', 'southwest')
title('\DeltaWS (20L-20F) - ERA 20C','Position',[0 0.7 1])
clear i; clear t; clear rr; clear tt; clear ans; clear p 
clear u2med;clear v2med
clear ws_reader; clear latr; clear lonr; clear wsmed
 %% comparar as reanalises (diferenca entre as medias)

wsmedtotal1 = mean(wsnoaa,3); wsmedtotal = mean(wsera,3);
u2medtotal = mean(u2,3); v2medtotal = mean(v2,3);
u1medtotal = mean(u1,3); v1medtotal = mean(v1,3);
clear u1; clear v1; clear u2; clear v2
clear wsnoaa; clear wsera

difrean = wsmedtotal - wsmedtotal1;

jjj = (-5:0.1:8);
figure('color',[1 1 1],'position',[10 305 940 815]);
m_proj('stereographic','lat',-90,'long',0,'radius',35);
colormap(cmocean('diff'))
m_contourf(lon,lat,difrean',jjj); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s')   
m_contour(lon,lat,difrean',jjj); hold on   
m_coast('color','k','linewidth',1.5); caxis([-7 7]); 
m_grid('xtick',[0 45 90 135 180 -45 -90 -135],'tickdir','out','ytick',...
    [-60 -75 -90],'tickdir','out', 'xaxisLocation', 'top', ...
    'yaxisLocation', 'middle','fontsize',7,'linestyle',':','linewidth',0.5);
title('ERA 20C vs NOAA 20CR (1900:2010)','Position',[0 0.7 1])

clear difrean; clear jjj

%% COASTAL -> MEAN 

figure('color',[1 1 1],'position',[10 305 940 515]); colormap(cmocean('dense')); 
subplot(121)
m_proj('lambert','lon',[0 60],'lat',[-80 -60]);
m_contourf(lon,lat,wsmedtotal1',ttt); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,wsmedtotal1',ttt);
m_coast('color','k','linewidth',1.5); 
m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),u1medtotal(1:2:end,1:2:end)',...
    v1medtotal(1:2:end,1:2:end)',2,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([0 12]); title('0 - 60E (NOAA 20CR)','Position',[0 0.23 1])
subplot(122)
m_proj('lambert','lon',[0 60],'lat',[-80 -60]);
m_contourf(lon,lat,wsmedtotal',ttt); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,wsmedtotal',ttt);
m_coast('color','k','linewidth',1.5); 
m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),u2medtotal(1:2:end,1:2:end)',...
    v2medtotal(1:2:end,1:2:end)',2,'color',[0.4 0.4 0.4]);

m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([0 12]); title('0 - 60E (ERA 20C)','Position',[0 0.23 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('color',[1 1 1],'position',[10 305 940 515]); colormap(cmocean('dense')); 
subplot(121)
m_proj('lambert','lon',[60 120],'lat',[-80 -60]);
m_contourf(lon,lat,wsmedtotal1',ttt); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,wsmedtotal1',ttt);
m_coast('color','k','linewidth',1.5); 
m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),u1medtotal(1:2:end,1:2:end)',...
    v1medtotal(1:2:end,1:2:end)',2,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([0 12]); title('60E - 120E (NOAA 20CR)','Position',[0 0.23 1])
subplot(122)
m_proj('lambert','lon',[60 120],'lat',[-80 -60]);
m_contourf(lon,lat,wsmedtotal',ttt); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,wsmedtotal',ttt);
m_coast('color','k','linewidth',1.5); 
m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),u2medtotal(1:2:end,1:2:end)',...
    v2medtotal(1:2:end,1:2:end)',2,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([0 12]); title('60E - 120E (ERA 20C)','Position',[0 0.23 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('color',[1 1 1],'position',[10 305 940 515]); colormap(cmocean('dense')); 
subplot(121)
m_proj('lambert','lon',[120 180],'lat',[-80 -60]);
m_contourf(lon,lat,wsmedtotal1',ttt); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,wsmedtotal1',ttt);
m_coast('color','k','linewidth',1.5); 
m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),u1medtotal(1:2:end,1:2:end)',...
    v1medtotal(1:2:end,1:2:end)',2,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([0 12]); title('120E - 180E (NOAA 20CR)','Position',[0 0.23 1])
subplot(122)
m_proj('lambert','lon',[120 180],'lat',[-80 -60]);
m_contourf(lon,lat,wsmedtotal',ttt); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,wsmedtotal',ttt);
m_coast('color','k','linewidth',1.5); 
m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),u2medtotal(1:2:end,1:2:end)',...
    v2medtotal(1:2:end,1:2:end)',2,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([0 12]); title('120E - 180E (ERA 20C)','Position',[0 0.23 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('color',[1 1 1],'position',[10 305 940 515]); colormap(cmocean('dense')); 
subplot(121)
m_proj('lambert','lon',[-180 -120],'lat',[-80 -60]);
m_contourf(lon,lat,wsmedtotal1',ttt); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,wsmedtotal1',ttt);
m_coast('color','k','linewidth',1.5); 
m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),u1medtotal(1:2:end,1:2:end)',...
    v1medtotal(1:2:end,1:2:end)',2,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([0 12]); title('180W - 120W (NOAA 20CR)','Position',[0 0.23 1])
subplot(122)
m_proj('lambert','lon',[-180 -120],'lat',[-80 -60]);
m_contourf(lon,lat,wsmedtotal',ttt); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,wsmedtotal',ttt);
m_coast('color','k','linewidth',1.5); 
m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),u2medtotal(1:2:end,1:2:end)',...
    v2medtotal(1:2:end,1:2:end)',2,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([0 12]); title('180W - 120W (ERA 20C)','Position',[0 0.23 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('color',[1 1 1],'position',[10 305 940 515]); colormap(cmocean('dense')); 
subplot(121)
m_proj('lambert','lon',[-120 -60],'lat',[-80 -60]);
m_contourf(lon,lat,wsmedtotal1',ttt); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,wsmedtotal1',ttt);
m_coast('color','k','linewidth',1.5); 
m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),u1medtotal(1:2:end,1:2:end)',...
    v1medtotal(1:2:end,1:2:end)',2,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([0 12]); title('120W - 60W (NOAA 20CR)','Position',[0 0.23 1])
subplot(122)
m_proj('lambert','lon',[-120 -60],'lat',[-80 -60]);
m_contourf(lon,lat,wsmedtotal',ttt); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,wsmedtotal',ttt);
m_coast('color','k','linewidth',1.5); 
m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),u2medtotal(1:2:end,1:2:end)',...
    v2medtotal(1:2:end,1:2:end)',2,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([0 12]); title('120W - 60W (ERA 20C)','Position',[0 0.23 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('color',[1 1 1],'position',[10 305 940 515]); colormap(cmocean('dense')); 
subplot(121)
m_proj('lambert','lon',[-60 0],'lat',[-80 -60]);
m_contourf(lon,lat,wsmedtotal1',ttt); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,wsmedtotal1',ttt);
m_coast('color','k','linewidth',1.5); 
m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),u1medtotal(1:2:end,1:2:end)',...
    v1medtotal(1:2:end,1:2:end)',2,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([0 12]); title('60W - 0 (NOAA 20CR)','Position',[0 0.23 1])
subplot(122)
m_proj('lambert','lon',[-60 0],'lat',[-80 -60]);
m_contourf(lon,lat,wsmedtotal',ttt); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,wsmedtotal',ttt);
m_coast('color','k','linewidth',1.5); 
m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),u2medtotal(1:2:end,1:2:end)',...
    v2medtotal(1:2:end,1:2:end)',2,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([0 12]); title('60W - 0 (ERA 20C)','Position',[0 0.23 1])

clear wsmedtotal; clear u2medtotal; clear v2medtotal; clear ans
clear wsmedtotal1; clear u1medtotal; clear v1medtotal; clear ans
clear ttt
%% COASTAL -> DIFF

figure('color',[1 1 1],'position',[10 305 940 515]); colormap(m_colmap('diverging'))
subplot(121)
m_proj('lambert','lon',[0 60],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsnoaa',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsnoaa',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u1(1:2:end,1:2:end)',...
    dif_v1(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 0 - 60E (NOAA 20CR)','Position',[0 0.23 1])
subplot(122)
m_proj('lambert','lon',[0 60],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsera',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsera',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u2(1:2:end,1:2:end)',...
    dif_v2(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 0 - 60E (ERA 20C)','Position',[0 0.23 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('color',[1 1 1],'position',[10 305 940 515]); colormap(m_colmap('diverging'))
subplot(121)
m_proj('lambert','lon',[60 120],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsnoaa',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsnoaa',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u1(1:2:end,1:2:end)',...
    dif_v1(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 60E - 120E (NOAA 20CR)','Position',[0 0.23 1])
subplot(122)
m_proj('lambert','lon',[60 120],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsera',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsera',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u2(1:2:end,1:2:end)',...
    dif_v2(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 60E - 120E (ERA 20C)','Position',[0 0.23 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('color',[1 1 1],'position',[10 305 940 515]); colormap(m_colmap('diverging'))
subplot(121)
m_proj('lambert','lon',[120 180],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsnoaa',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsnoaa',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u1(1:2:end,1:2:end)',...
    dif_v1(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 120E - 180E (NOAA 20CR)','Position',[0 0.23 1])
subplot(122)
m_proj('lambert','lon',[120 180],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsera',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsera',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u2(1:2:end,1:2:end)',...
    dif_v2(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 120E - 180E (ERA 20C)','Position',[0 0.23 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('color',[1 1 1],'position',[10 305 940 515]); colormap(m_colmap('diverging'))
subplot(121)
m_proj('lambert','lon',[-180 -120],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsnoaa',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsnoaa',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u1(1:2:end,1:2:end)',...
    dif_v1(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 180W - 120W (NOAA 20CR)','Position',[0 0.23 1])
subplot(122)
m_proj('lambert','lon',[-180 -120],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsera',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsera',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u2(1:2:end,1:2:end)',...
    dif_v2(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 180W - 120W (ERA 20C)','Position',[0 0.23 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('color',[1 1 1],'position',[10 305 940 515]); colormap(m_colmap('diverging'))
subplot(121)
m_proj('lambert','lon',[-120 -60],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsnoaa',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsnoaa',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u1(1:2:end,1:2:end)',...
    dif_v1(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 120W - 60W (NOAA 20CR)','Position',[0 0.23 1])
subplot(122)
m_proj('lambert','lon',[-120 -60],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsera',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsera',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u2(1:2:end,1:2:end)',...
    dif_v2(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 120W - 60W (ERA 20C)','Position',[0 0.23 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('color',[1 1 1],'position',[10 305 940 515]); colormap(m_colmap('diverging'))
subplot(121)
m_proj('lambert','lon',[-60 0],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsnoaa',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsnoaa',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u1(1:2:end,1:2:end)',...
    dif_v1(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 60W - 0 (NOAA 20CR)','Position',[0 0.23 1])
subplot(122)
m_proj('lambert','lon',[-60 0],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsera',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsera',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u2(1:2:end,1:2:end)',...
    dif_v2(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 60W - 0 (ERA 20C)','Position',[0 0.23 1])


clear dif_wsnoaa; clear dif_u1; clear dif_v1; 
clear dif_wsera; clear dif_u2; clear dif_v2;
clear ans; clear x; clear y; clear p; clear pp; clear t; clear tt


%%%%%%%%%%%%
figure('color',[1 1 1],'position',[10 305 940 515]); colormap(m_colmap('diverging'))
subplot(121)
m_proj('lambert','lon',[90 135],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsnoaa',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsnoaa',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u1(1:2:end,1:2:end)',...
    dif_v1(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 90-135^oE (NOAA 20CR)','Position',[0 0.23 1])
subplot(122)
m_proj('lambert','lon',[90 135],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsera',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsera',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u2(1:2:end,1:2:end)',...
    dif_v2(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 90-135^oE (ERA 20C)','Position',[0 0.23 1])


figure('color',[1 1 1],'position',[10 305 940 515]); colormap(m_colmap('diverging'))
subplot(121)
m_proj('lambert','lon',[-135 -90],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsnoaa',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsnoaa',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u1(1:2:end,1:2:end)',...
    dif_v1(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 135-90^oW (NOAA 20CR)','Position',[0 0.23 1])
subplot(122)
m_proj('lambert','lon',[-135 -90],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsera',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsera',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u2(1:2:end,1:2:end)',...
    dif_v2(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 135-90^oW (ERA 20C)','Position',[0 0.23 1])


figure('color',[1 1 1],'position',[10 305 940 515]); colormap(m_colmap('diverging'))
subplot(121)
m_proj('lambert','lon',[-90 -45],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsnoaa',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsnoaa',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u1(1:2:end,1:2:end)',...
    dif_v1(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 90-45^oW (NOAA 20CR)','Position',[0 0.23 1])
subplot(122)
m_proj('lambert','lon',[-90 -45],'lat',[-80 -60]);
m_contourf(lon,lat,dif_wsera',pp); hold on
t = colorbar('Location', 'eastoutside'); tt = get(t,'title'); set(tt,'string','m/s') 
m_contour(lon,lat,dif_wsera',pp)
m_coast('color','k','linewidth',1.5); 
p = m_quiver(x(1:2:end,1:2:end),y(1:2:end,1:2:end),dif_u2(1:2:end,1:2:end)',...
    dif_v2(1:2:end,1:2:end)',1,'color',[0.5 0.5 0.5]);
m_grid('xtick',4,'ytick',4,'tickdir','out','box','fancy','tickdir','in','xaxisLocation', 'top');
caxis([-2.5 2.5]); title('\DeltaWS 90-45^oW (ERA 20C)','Position',[0 0.23 1])












%% MAPAS SETORES
a = linspace(-180,180,7);
bndry_lat =[-60*ones(1,61) -80*ones(1,61)]; % Outline of a box

for j = 1

    bndry_lon=[a(j):1:a(j+1) a(j+1):-1:a(j)];

    figure('color',[1 1 1],'position',[108 305 830 817]);
    m_proj('stereographic','lat',-90,'long',0,'radius',40);
    m_coast('color',[0.5 0.5 0.5],'linewidth',1.5);
    m_grid('ytick',[-60 -68 -80 -90],'xtick',[-120 -60],...
        'tickdir','out', 'xaxisLocation', 'top'...
        ,'box','on', 'fontsize',12,'linewidth',1.5);
    m_hatch(bndry_lon, bndry_lat,'single',30,5,'color',[0.3020 0.7451 0.9333],'linewidth', 2);

end














