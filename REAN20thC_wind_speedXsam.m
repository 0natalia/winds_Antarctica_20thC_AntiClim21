%__________________________________________________________________________
%                       AntClim21 Workshop - BAS                          %
%               Wind changes in Antarctica during the 20th C              %
%              Reanalysis data from NOAA-20CR e ECMWF-ERA20C              %

% This script reads wind and slp data; calculates the SAM index and its
% trends for before and after 1975; calculates and plots the correlation
% beteen SAM index and winds

% Nat?lia Silva; natalia3.silva@usp.br
% 2020
%__________________________________________________________________________

clear all; close all

%% DADOS VENTO

%%% ERA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
era20u = ('home/natalia/reanalises/ERA20C/era_uwnd.nc');
era20v = ('home/natalia/reanalises/ERA20C/era_vwnd.nc');
lat = ncread(era20u,'lat'); l = lat<-54;
lon = ncread(era20u,'lon'); 
time2 = ncread(era20u,'time'); tempo = datenum(datetime(1900,1,1)+hours(time2)); clear time2
meses = linspace(1,length(tempo),length(tempo));

u2 = ncread(era20u,'u'); v2 = ncread(era20v,'v'); % nivel 5 = superficie
u2 = squeeze(u2(:,l,5,:)); v2 = squeeze(v2(:,l,5,:));
wsera = hypot(u2,v2);

wsera = cat(1,wsera,wsera(end,:,:)); lon(end+1) = 180;
wsera = [wsera(:,1,:) wsera]; lat = [-90; lat];

[wsera_anual,tanos] = downsample_ts(wsera,tempo,'year');

%%% NOAA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noaav = ('home/natalia/reanalises/NOAA20CR/vwnd.10m.mon.mean.nc');
noaau = ('home/natalia/reanalises/NOAA20CR/uwnd.10m.mon.mean.nc');
v1 = ncread(noaav,'vwnd'); v1 = v1(:,l,349:1680);
u1 = ncread(noaau,'uwnd'); u1 = u1(:,l,349:1680);
wsnoaa = hypot(u1,v1);

wsnoaa = cat(1,wsnoaa,wsnoaa(end,:,:));
wsnoaa = [wsnoaa(:,1,:) wsnoaa];

[wsnoaa_anual,~] = downsample_ts(wsnoaa,tempo,'year');

clear era20u; clear era20v; clear noaau; clear noaav
clear v1; clear u1; clear v2; clear u2; clear l

%% DADOS PRESSAO

%%% ERA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pera = ('home/natalia/reanalises/ERA20C/era_slp.nc');
slpera = ncread(pera,'msl'); slpera = slpera/100; % hPa
slpera = cat(1,slpera, slpera(end,:,:));
slpera = [slpera(:,1,:) slpera];

%%% NOAA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pnoaa = ('home/natalia/reanalises/NOAA20CR/slp.mon.mean.nc');
slpnoaa = ncread(pnoaa,'prmsl'); 
slpnoaa = slpnoaa(:,:,589:1920); slpnoaa = slpnoaa/100; % hPa
slpnoaa = cat(1, slpnoaa, slpnoaa(end,:,:));
slpnoaa = [slpnoaa(:,1,:) slpnoaa];

clear pera; clear pnoaa

 %% SAM ANUAL

n_bndry = lat>-42 & lat<-39; s_bndry = lat<-64 & lat> -66.5;
slpnoaa40a = slpnoaa(:,n_bndry,:); slpnoaa65a = slpnoaa(:,s_bndry,:);
slpera40a = slpera(:,n_bndry,:); slpera65a = slpera(:,s_bndry,:);

[slpnoaa40a,~] = downsample_ts(slpnoaa40a,tempo,'year');
[slpnoaa65a,~] = downsample_ts(slpnoaa65a,tempo,'year');
[slpera40a,~] = downsample_ts(slpera40a,tempo,'year');
[slpera65a,~] = downsample_ts(slpera65a,tempo,'year');

slpnoaa40a = squeeze(mean(slpnoaa40a,2)); slpnoaa65a = squeeze(mean(slpnoaa65a,2));
slpera40a = squeeze(mean(slpera40a,2)); slpera65a = squeeze(mean(slpera65a,2));

% normalizar
nor_noaa40a = (squeeze(mean(slpnoaa40a)) - mean2(slpnoaa40a))/std(mean(slpnoaa40a));
nor_noaa65a = (squeeze(mean(slpnoaa65a)) - mean2(slpnoaa65a))/std(mean(slpnoaa65a));
nor_era40a = (squeeze(mean(slpera40a)) - mean2(slpera40a))/std(mean(slpera40a));
nor_era65a = (squeeze(mean(slpera65a)) - mean2(slpera65a))/std(mean(slpera65a));

sam_noaa_anual = nor_noaa40a - nor_noaa65a; 
sam_era_anual = nor_era40a - nor_era65a;

clear nor_noaa40a; clear nor_noaa65a; clear nor_era40a; clear nor_era65a
clear slpnoaa40a; clear slpnoaa65a; clear slpera40a; clear slpera65a


%% SAM MENSAL

slpnoaa40m = slpnoaa(:,n_bndry,:); slpnoaa65m = slpnoaa(:,s_bndry,:);
slpera40m = slpera(:,n_bndry,:); slpera65m = slpera(:,s_bndry,:);

clear n_bndry; clear s_bndry;

slpnoaa40m = squeeze(mean(slpnoaa40m,2)); slpnoaa65m = squeeze(mean(slpnoaa65m,2));
slpera40m = squeeze(mean(slpera40m,2)); slpera65m = squeeze(mean(slpera65m,2));

% normalizar
nor_noaa40m = (squeeze(mean(slpnoaa40m)) - mean2(slpnoaa40m))/std(mean(slpnoaa40m));
nor_noaa65m = (squeeze(mean(slpnoaa65m)) - mean2(slpnoaa65m))/std(mean(slpnoaa65m));
nor_era40m = (squeeze(mean(slpera40m)) - mean2(slpera40m))/std(mean(slpera40m));
nor_era65m = (squeeze(mean(slpera65m)) - mean2(slpera65m))/std(mean(slpera65m));

sam_noaa = nor_noaa40m - nor_noaa65m; 
sam_era = nor_era40m - nor_era65m;

clear nor_noaa40m; clear nor_noaa65m; clear nor_era40m; clear nor_era65m
clear slpnoaa40m; clear slpnoaa65m; clear slpera40m; clear slpera65m
%clear slpera; clear slpnoaa

%% SAM TREND
lsq_noaa = fitlm(meses,sam_noaa);
yo_noaa = lsq_noaa.Coefficients.Estimate(1);
trend_noaa = lsq_noaa.Coefficients.Estimate(2); erro_noaa = lsq_noaa.Coefficients.SE(2);
line_noaa = yo_noaa+trend_noaa*meses; trend_noaa = (trend_noaa*120); erro_noaa = erro_noaa*120;
p_noaa = lsq_noaa.Coefficients.pValue(2); 

if p_noaa < 0.01
    SIG_noaa = ('p < 0.01');
elseif p_noaa < 0.05
    SIG_noaa = ('p < 0.05');
elseif p_noaa < 0.1
    SIG_noaa = ('p < 0.1');
else
    SIG_noaa = ('p > 0.1');
end
strnoaa = mat2str(trend_noaa); strErnoaa = mat2str(erro_noaa); 
disp(['NOAA Tend: ',strnoaa,' \pm ',strErnoaa]); disp(SIG_noaa)

clear strnoaa; clear strErnoaa; clear SIG_noaa; clear p_noaa; 
clear trend_noaa; clear erro_noaa; clear lsq_noaa; clear yo_noaa;

lsq_era = fitlm(meses,sam_era);
yo_era = lsq_era.Coefficients.Estimate(1);
trend_era = lsq_era.Coefficients.Estimate(2); erro_era = lsq_era.Coefficients.SE(2);
line_era = yo_era+trend_era*meses; trend_era = (trend_era*120); erro_era = erro_era*120;
p_era = lsq_era.Coefficients.pValue(2); 

if p_era < 0.01
    SIG_era = ('p < 0.01');
elseif p_era < 0.05
    SIG_era = ('p < 0.05');
elseif p_era < 0.1
    SIG_era = ('p < 0.1');
else
    SIG_era = ('p > 0.1');
end
strera = mat2str(trend_era); strErera = mat2str(erro_era);
disp(['Eram Trend: ',strera,' \pm ',strErera]); disp(SIG_era)

clear strera; clear strErera; clear SIG_era; clear p_era; 
clear trend_era; clear erro_era; clear lsq_era; clear yo_era;

% %******** SAM NOAA B75 **************************

lsq_noaa_B75 = fitlm(meses(1:901),sam_noaa(1:901));
yo_noaa_B75 = lsq_noaa_B75.Coefficients.Estimate(1);
trend_noaa_B75 = lsq_noaa_B75.Coefficients.Estimate(2); 
erro_noaa_B75 = lsq_noaa_B75.Coefficients.SE(2);
line_noaa_B75 = yo_noaa_B75+trend_noaa_B75*meses(1:901); trend_noaa_B75 = (trend_noaa_B75*120); 
erro_noaa_B75 = erro_noaa_B75*120;
p_noaa_B75 = lsq_noaa_B75.Coefficients.pValue(2); 

if p_noaa_B75 < 0.01
    SIG_noaa_B75 = ('p < 0.01');
elseif p_noaa_B75 < 0.05
    SIG_noaa_B75 = ('p < 0.05');
elseif p_noaa_B75 < 0.1
    SIG_noaa_B75 = ('p < 0.1');
else
    SIG_noaa_B75 = ('p > 0.1');
end
strnoaa_B75 = mat2str(trend_noaa_B75); strErnoaa_B75 = mat2str(erro_noaa_B75); 
disp(['NOAA_B75 Tend: ',strnoaa_B75,' \pm ',strErnoaa_B75]); disp(SIG_noaa_B75)

clear strnoaa_B75; clear strErnoaa_B75; clear SIG_noaa_B75; clear p_noaa_B75; 
clear trend_noaa_B75; clear erro_noaa_B75; clear lsq_noaa_B75; clear yo_noaa_B75;

lsq_era_B75 = fitlm(meses(1:901),sam_era(1:901));
yo_era_B75 = lsq_era_B75.Coefficients.Estimate(1);
trend_era_B75 = lsq_era_B75.Coefficients.Estimate(2); erro_era_B75 = lsq_era_B75.Coefficients.SE(2);
line_era_B75 = yo_era_B75+trend_era_B75*meses(1:901); trend_era_B75 = (trend_era_B75*120); 
erro_era_B75 = erro_era_B75*120;
p_era_B75 = lsq_era_B75.Coefficients.pValue(2); 

if p_era_B75 < 0.01
    SIG_era_B75 = ('p < 0.01');
elseif p_era_B75 < 0.05
    SIG_era_B75 = ('p < 0.05');
elseif p_era_B75 < 0.1
    SIG_era_B75 = ('p < 0.1');
else
    SIG_era_B75 = ('p > 0.1');
end
strera_B75 = mat2str(trend_era_B75); strErera_B75 = mat2str(erro_era_B75);
disp(['Eram_B75 Trend: ',strera_B75,' \pm ',strErera_B75]); disp(SIG_era_B75)

clear strera_B75; clear strErera_B75; clear SIG_era_B75; clear p_era_B75; 
clear trend_era_B75; clear erro_era_B75; clear lsq_era_B75; clear yo_era_B75;

%******** SAM NOAA A75 **************************

lsq_noaa_A75 = fitlm(meses(901:end),sam_noaa(901:end));
yo_noaa_a75 = lsq_noaa_A75.Coefficients.Estimate(1);
trend_noaa_A75 = lsq_noaa_A75.Coefficients.Estimate(2); 
erro_noaa_A75 = lsq_noaa_A75.Coefficients.SE(2);
line_noaa_A75 = yo_noaa_a75+trend_noaa_A75*meses(901:end); trend_noaa_A75 = (trend_noaa_A75*120); 
erro_noaa_A75 = erro_noaa_A75*120;
p_noaa_A75 = lsq_noaa_A75.Coefficients.pValue(2); 

if p_noaa_A75 < 0.01
    SIG_noaa_A75 = ('p < 0.01');
elseif p_noaa_A75 < 0.05
    SIG_noaa_75 = ('p < 0.05');
elseif p_noaa_A75 < 0.1
    SIG_noaa_75 = ('p < 0.1');
else
    SIG_noaa_75 = ('p > 0.1');
end
strnoaa_A75 = mat2str(trend_noaa_A75); strErnoaa_A75 = mat2str(erro_noaa_A75); 
disp(['NOAA_A75 Tend: ',strnoaa_A75,' \pm ',strErnoaa_A75]); disp(SIG_noaa_A75)

clear strnoaa_A75; clear strErnoaa_A75; clear SIG_noaa_A75; clear p_noaa_A75; 
clear trend_noaa_A75; clear erro_noaa_A75; clear lsq_noaa_A75; clear yo_noaa_A75;

lsq_era_A75 = fitlm(meses(901:end),sam_era(901:end));
yo_era_A75 = lsq_era_A75.Coefficients.Estimate(1);
trend_era_A75 = lsq_era_A75.Coefficients.Estimate(2); erro_era_A75 = lsq_era_A75.Coefficients.SE(2);
line_era_A75 = yo_era_A75+trend_era_A75*meses(901:end); trend_era_A75 = (trend_era_A75*120); 
erro_era_A75 = erro_era_A75*120;
p_era_A75 = lsq_era_A75.Coefficients.pValue(2); 

if p_era_A75 < 0.01
    SIG_era_A75 = ('p < 0.01');
elseif p_era_A75 < 0.05
    SIG_era_A75 = ('p < 0.05');
elseif p_era_A75 < 0.1
    SIG_era_A75 = ('p < 0.1');
else
    SIG_era_A75 = ('p > 0.1');
end
strera_A75 = mat2str(trend_era_A75); strErera_A75 = mat2str(erro_era_A75);
disp(['Eram_B75 Trend: ',strera_A75,' \pm ',strErera_A75]); disp(SIG_era_A75)

clear strera_A75; clear strErera_A75; clear SIG_era_A75; clear p_era_A75; 
clear trend_era_A75; clear erro_era_A75; clear lsq_era_A75; clear yo_era_A75;

%%%%%%%%%%%%%%%%%%%%%%% PLOT
figure('color',[1 1 1],'position',[108 305 900 500]) 
subplot(2,1,1); plot(tempo,sam_noaa,'color', [0.87 0.87 0.87],'linewidth',1);
hold on; plot(tanos,sam_noaa_anual,'k','linewidth',1.3);
title('SAM (NOAA20CR)'); xlabel('Time');
datetick('x',10,'keepticks'); xlim([tanos(1) tanos(end)]); ylim([-5 5])
plot([tempo(1) tempo(end)],[0 0],'Color',[0.6 0.6 0.6],'LineStyle','--')
text(tanos(2),4,'Trend_{B75}: 0.054 \pm 0.022; p < 0.05')
text(tanos(2),3,'Trend_{A75}: 0.350 \pm 0.078; p < 0.01')
plot(tempo(1:901),line_noaa_B75,'r','linewidth',1.7) 
plot(tempo(901:end),line_noaa_A75,'r','linewidth',1.7) 
plot([tanos(75) tanos(75)], [-7 7],'LineStyle','--', 'color', [0.4 0.4 0.4])


subplot(2,1,2); plot(tempo,sam_era,'color',[0.87 0.87 0.87], 'linewidth',1);
hold on; plot(tanos,sam_era_anual,'k','linewidth',1.3);
title('SAM (ERA20C)'); xlabel('Time');
datetick('x',10,'keepticks'); xlim([tanos(1) tanos(end)]); ylim([-5 5])
plot([tempo(1) tempo(end)],[0 0],'Color',[0.6 0.6 0.6],'LineStyle','--')
text(tanos(2),4,'Trend_{B75}: 0.276 \pm 0.024; p < 0.01')
text(tanos(2),3,'Trend_{A75}: 0.408 \pm 0.063; p < 0.01')
plot(tempo(1:901),line_era_B75,'r','linewidth',1.7)
plot(tempo(901:end),line_era_A75,'r','linewidth',1.7)
plot([tanos(75) tanos(75)], [-7 7],'LineStyle','--', 'color', [0.4 0.4 0.4])

clear line_noaa; clear line_era

%% CORRELACAO
latc = lat(lat<-54);
corr_e = zeros(193,20); p_e = zeros(193,20);
corr_n = zeros(193,20); p_n = zeros(193,20);

% ERA
for i = 1:length(lon)
    for j = 1:length(latc)
        [ce, pe] = corrcoef(sam_era, wsera(i,j,:),'alpha',0.01);
        [cn, pn] = corrcoef(sam_noaa, wsnoaa(i,j,:),'alpha',0.01);
        
        corr_e(i,j) = ce(2,1); p_e(i,j) = pe(2,1);
        corr_n(i,j) = cn(2,1); p_n(i,j) = pn(2,1);
    end
end
clear ce; clear pe; clear cn; clear pn
clear i; clear j

% corr significativa ao nivel de 1%
pp = p_e < 0.01; p_e(~pp) = NaN; ppp = p_n < 0.01; p_n(~ppp) = NaN;
clear pp; clear ppp

% CORR ERA 20C
v=(-1:0.01:1);
figure('color',[1 1 1],'position',[10 805 800 700]); 
colormap(flipud(cbrewer('div','Spectral',40)))
m_proj('stereographic','lat',-90,'long',0,'radius',35);
m_contourf(lon,latc,corr_e',v); hold on
m_contour(lon,latc,corr_e',v)
m_contour(lon,latc,p_e',':k')
%[c] = m_contour(lon,lat,mean(slpera,3)',[980, 990,1000,1010, 1020], 'color',[0.7 0.7 0.7]);
m_grid('xtick',[0 45 90 135 180 -45 -90 -135],'tickdir','out','ytick',...
    [-60 -75 -90],'tickdir','out', 'xaxisLocation', 'top', ...
    'yaxisLocation', 'middle','fontsize',7,'linestyle',':','linewidth',0.5);
m_coast('color','k','linewidth',1.5); caxis([-1 1])
title('Corr SAM x WS (ERA 20C)','Position',[0 0.7 1]);
t = colorbar; tt = get(t,'title'); set(tt,'string','r')
clear t; clear ttt; clear tt; clear ans;
%clabel(c, 'manual','color', [0.5 0.5 0.5])

% CORR ERA 20C
v=(-1:0.01:1);
figure('color',[1 1 1],'position',[10 805 800 700]); 
colormap(flipud(cbrewer('div','Spectral',40)))
m_proj('stereographic','lat',-90,'long',0,'radius',35);
m_contourf(lon,latc,corr_n',v); hold on
m_contour(lon,latc,corr_n',v)
m_contour(lon,latc,p_n',':k')
%[cc] = m_contour(lon,lat,mean(slpnoaa,3)',[980, 990,1000,1010, 1015, 1020], 'color',[0.7 0.7 0.7]);
m_grid('xtick',[0 45 90 135 180 -45 -90 -135],'tickdir','out','ytick',...
    [-60 -75 -90],'tickdir','out', 'xaxisLocation', 'top', ...
    'yaxisLocation', 'middle','fontsize',7,'linestyle',':','linewidth',0.5);
m_coast('color','k','linewidth',1.5); caxis([-1 1])
title('Corr SAM x WS (NOAA 20CR)','Position',[0 0.7 1]);
t = colorbar; tt = get(t,'title'); set(tt,'string','r')
clear t; clear ttt; clear tt; clear ans; clear v
%clabel(cc, 'manual')


%% MENSAL 58to00
% tempo_5800 = tempo(697:1212); meses_5800 = meses(697:1212);
% n_bndry = lat>-42 & lat<-39; s_bndry = lat<-64 & lat> -66.5;
% 
% slpnoaa40m_5800 = slpnoaa(:,n_bndry,697:1212); slpnoaa65m_5800 = slpnoaa(:,s_bndry,697:1212);
% slpera40m_5800 = slpera(:,n_bndry,697:1212); slpera65m_5800 = slpera(:,s_bndry,697:1212);
% 
% slpnoaa40m_5800 = squeeze(mean(slpnoaa40m_5800,2)); slpnoaa65m_5800 = squeeze(mean(slpnoaa65m_5800,2));
% slpera40m_5800 = squeeze(mean(slpera40m_5800,2)); slpera65m_5800 = squeeze(mean(slpera65m_5800,2));
% 
% % normalizar
% nor_noaa40m_5800 = (squeeze(mean(slpnoaa40m_5800)) - mean2(slpnoaa40m_5800))/std(mean(slpnoaa40m_5800));
% nor_noaa65m_5800 = (squeeze(mean(slpnoaa65m_5800)) - mean2(slpnoaa65m_5800))/std(mean(slpnoaa65m_5800));
% nor_era40m_5800 = (squeeze(mean(slpera40m_5800)) - mean2(slpera40m_5800))/std(mean(slpera40m_5800));
% nor_era65m_5800 = (squeeze(mean(slpera65m_5800)) - mean2(slpera65m_5800))/std(mean(slpera65m_5800));
% 
% sam_noaa_5800 = nor_noaa40m_5800 - nor_noaa65m_5800; 
% sam_era_5800 = nor_era40m_5800 - nor_era65m_5800;
% 
% clear nor_noaa40m_5800; clear nor_noaa65m_5800; clear nor_era40m_5800; clear nor_era65m_5800
% clear slpnoaa40m_5800; clear slpnoaa65m_5800; clear slpera40m_5800; clear slpera65m_5800
% 
% lsqnoaam_5800 = fitlm(meses_5800,sam_noaa_5800);
% yonoaam_5800 = lsqnoaam_5800.Coefficients.Estimate(1);
% trendnoaam_5800 = lsqnoaam_5800.Coefficients.Estimate(2);
% erronoaam_5800 = lsqnoaam_5800.Coefficients.SE(2);
% linenoaam_5800 = yonoaam_5800+trendnoaam_5800*meses_5800; 
% trendnoaam_5800 = (trendnoaam_5800*120); erronoaam_5800 = erronoaam_5800*120;
% p_valuenoaam_5800 = lsqnoaam_5800.Coefficients.pValue(2); 
% 
% if p_valuenoaam_5800 < 0.01
%     SIGnoaam_5800 = ('p < 0.01');
% elseif p_valuenoaam_5800 < 0.05
%     SIGnoaam_5800 = ('p < 0.05');
% elseif p_valuenoaam_5800 < 0.1
%     SIGnoaam_5800 = ('p < 0.1');
% else
%     SIGnoaam_5800 = ('p > 0.1');
% end
% strnoaam = mat2str(trendnoaam_5800); strErnoaam = mat2str(erronoaam_5800); 
% disp(['NOAA 5800 Tend: ',strnoaam,' \pm ',strErnoaam,'(ms^-^1/decade)']); disp(SIGnoaam_5800)
% clear strnoaam; clear strErnoaam; clear SIGnoaam_5800; clear p_valuenoaam_5800; 
% clear trendnoaam_5800; clear erronoaam_5800; clear lsqnoaam_5800; clear yonoaamy_5800;
% 
% lsqeram_5800 = fitlm(meses_5800,sam_era_5800);
% yoeram_5800 = lsqeram_5800.Coefficients.Estimate(1);
% trenderam_5800 = lsqeram_5800.Coefficients.Estimate(2);
% erroeram_5800 = lsqeram_5800.Coefficients.SE(2);
% lineeram_5800 = yoeram_5800+trenderam_5800*meses_5800; 
% trenderam_5800 = (trenderam_5800*120); erroeram_5800 = erroeram_5800*120;
% p_valueeram_5800 = lsqeram_5800.Coefficients.pValue(2); 
% 
% if p_valueeram_5800 < 0.01
%     SIGeram_5800 = ('p < 0.01');
% elseif p_valueeram_5800 < 0.05
%     SIGeram_5800 = ('p < 0.05');
% elseif p_valueeram_5800 < 0.1
%     SIGeram_5800 = ('p < 0.1');
% else
%     SIGeram_5800 = ('p > 0.1');
% end
% streram = mat2str(trenderam_5800); strEreram = mat2str(erroeram_5800); 
% disp(['ERA 5800 Tend: ',streram,' \pm ',strEreram]); disp(SIGeram_5800)
% clear streram; clear strEreram; clear SIGeram_5800; clear p_valueeram_5800; 
% clear trenderam_5800; clear erroeram_5800; clear lsqeram_5800; clear yoeram_5800;

%% MENSAL 82 - 98

% tempo_8298 = tempo(985:1188); meses_8298 = meses(985:1188);
% n_bndry = lat>-42 & lat<-39; s_bndry = lat<-64 & lat> -66.5;
% 
% slpnoaa40m_8298 = slpnoaa(:,n_bndry,985:1188); slpnoaa65m_8298 = slpnoaa(:,s_bndry,985:1188);
% slpera40m_8298 = slpera(:,n_bndry,985:1188); slpera65m_8298 = slpera(:,s_bndry,985:1188);
% 
% slpnoaa40m_8298 = squeeze(mean(slpnoaa40m_8298,2)); slpnoaa65m_8298 = squeeze(mean(slpnoaa65m_8298,2));
% slpera40m_8298 = squeeze(mean(slpera40m_8298,2)); slpera65m_8298 = squeeze(mean(slpera65m_8298,2));
% 
% % normalizar
% nor_noaa40m_8298 = (squeeze(mean(slpnoaa40m_8298)) - mean2(slpnoaa40m_8298))/std(mean(slpnoaa40m_8298));
% nor_noaa65m_8298 = (squeeze(mean(slpnoaa65m_8298)) - mean2(slpnoaa65m_8298))/std(mean(slpnoaa65m_8298));
% nor_era40m_8298 = (squeeze(mean(slpera40m_8298)) - mean2(slpera40m_8298))/std(mean(slpera40m_8298));
% nor_era65m_8298 = (squeeze(mean(slpera65m_8298)) - mean2(slpera65m_8298))/std(mean(slpera65m_8298));
% 
% sam_noaa_8298 = nor_noaa40m_8298 - nor_noaa65m_8298; 
% sam_era_8298 = nor_era40m_8298 - nor_era65m_8298;
% 
% clear nor_noaa40m_8298; clear nor_noaa65m_8298; clear nor_era40m_8298; clear nor_era65m_8298
% clear slpnoaa40m_8298; clear slpnoaa65m_8298; clear slpera40m_8298; clear slpera65m_8298
% 
% lsqnoaam_8298 = fitlm(meses_8298,sam_noaa_8298);
% yonoaam_8298 = lsqnoaam_8298.Coefficients.Estimate(1);
% trendnoaam_8298 = lsqnoaam_8298.Coefficients.Estimate(2);
% erronoaam_8298 = lsqnoaam_8298.Coefficients.SE(2);
% linenoaam_8298 = yonoaam_8298+trendnoaam_8298*meses_8298; 
% trendnoaam_8298 = (trendnoaam_8298*120); erronoaam_8298 = erronoaam_8298*120;
% p_valuenoaam_8298 = lsqnoaam_8298.Coefficients.pValue(2); 
% 
% if p_valuenoaam_8298 < 0.01
%     SIGnoaam_8298 = ('p < 0.01');
% elseif p_valuenoaam_8298 < 0.05
%     SIGnoaam_8298 = ('p < 0.05');
% elseif p_valuenoaam_8298 < 0.1
%     SIGnoaam_8298 = ('p < 0.1');
% else
%     SIGnoaam_8298 = ('p > 0.1');
% end
% strnoaam = mat2str(trendnoaam_8298); strErnoaam = mat2str(erronoaam_8298); 
% disp(['NOAA 8298 Tend: ',strnoaam,' \pm ',strErnoaam,'(ms^-^1/decade)']); disp(SIGnoaam_8298)
% clear strnoaam; clear strErnoaam; clear SIGnoaam_8298; clear p_valuenoaam_8298; 
% clear trendnoaam_8298; clear erronoaam_8298; clear lsqnoaam_8298; clear yonoaamy_8298;
% 
% lsqeram_8298 = fitlm(meses_8298,sam_era_8298);
% yoeram_8298 = lsqeram_8298.Coefficients.Estimate(1);
% trenderam_8298 = lsqeram_8298.Coefficients.Estimate(2);
% erroeram_8298 = lsqeram_8298.Coefficients.SE(2);
% lineeram_8298 = yoeram_8298+trenderam_8298*meses_8298; 
% trenderam_8298 = (trenderam_8298*120); erroeram_8298 = erroeram_8298*120;
% p_valueeram_8298 = lsqeram_8298.Coefficients.pValue(2); 
% 
% if p_valueeram_8298 < 0.01
%     SIGeram_8298 = ('p < 0.01');
% elseif p_valueeram_8298 < 0.05
%     SIGeram_8298 = ('p < 0.05');
% elseif p_valueeram_8298 < 0.1
%     SIGeram_8298 = ('p < 0.1');
% else
%     SIGeram_8298 = ('p > 0.1');
% end
% streram = mat2str(trenderam_8298); strEreram = mat2str(erroeram_8298); 
% disp(['ERA 8298 Tend: ',streram,' \pm ',strEreram]); disp(SIGeram_8298)
% clear streram; clear strEreram; clear SIGeram_8298; clear p_valueeram_8298; 
% clear trenderam_8298; clear erroeram_8298; clear lsqeram_8298; clear yoeram_8298;



%% MENSAL 65 - 97

% tempo_6597 = tempo(781:1176); meses_6597 = meses(781:1176);
% n_bndry = lat>-42 & lat<-39; s_bndry = lat<-64 & lat> -66.5;
% 
% slpnoaa40m_6597 = slpnoaa(:,n_bndry,781:1176); slpnoaa65m_6597 = slpnoaa(:,s_bndry,781:1176);
% slpera40m_6597 = slpera(:,n_bndry,781:1176); slpera65m_6597 = slpera(:,s_bndry,781:1176);
% 
% slpnoaa40m_6597 = squeeze(mean(slpnoaa40m_6597,2)); slpnoaa65m_6597 = squeeze(mean(slpnoaa65m_6597,2));
% slpera40m_6597 = squeeze(mean(slpera40m_6597,2)); slpera65m_6597 = squeeze(mean(slpera65m_6597,2));
% 
% % normalizar
% nor_noaa40m_6597 = (squeeze(mean(slpnoaa40m_6597)) - mean2(slpnoaa40m_6597))/std(mean(slpnoaa40m_6597));
% nor_noaa65m_6597 = (squeeze(mean(slpnoaa65m_6597)) - mean2(slpnoaa65m_6597))/std(mean(slpnoaa65m_6597));
% nor_era40m_6597 = (squeeze(mean(slpera40m_6597)) - mean2(slpera40m_6597))/std(mean(slpera40m_6597));
% nor_era65m_6597 = (squeeze(mean(slpera65m_6597)) - mean2(slpera65m_6597))/std(mean(slpera65m_6597));
% 
% sam_noaa_6597 = nor_noaa40m_6597 - nor_noaa65m_6597; 
% sam_era_6597 = nor_era40m_6597 - nor_era65m_6597;
% 
% clear nor_noaa40_6597; clear nor_noaa65_6597; clear nor_era40_6597; clear nor_era65_6597
% clear slpnoaa40_6597; clear slpnoaa65_6597; clear slpera40_6597; clear slpera65_6597
% 
% lsqnoaam_6597 = fitlm(meses_6597,sam_noaa_6597);
% yonoaa_6597 = lsqnoaam_6597.Coefficients.Estimate(1);
% trendnoaam_6597 = lsqnoaam_6597.Coefficients.Estimate(2);
% erronoaam_6597 = lsqnoaam_6597.Coefficients.SE(2);
% linenoaam_6597 = yonoaa_6597+trendnoaam_6597*meses_6597; 
% trendnoaam_6597 = (trendnoaam_6597*120); erronoaam_6597 = erronoaam_6597*120;
% p_valuenoaam_6597 = lsqnoaam_6597.Coefficients.pValue(2); 
% 
% if p_valuenoaam_6597 < 0.01
%     SIGnoaam_6597 = ('p < 0.01');
% elseif p_valuenoaam_6597 < 0.05
%     SIGnoaam_6597 = ('p < 0.05');
% elseif p_valuenoaam_6597 < 0.1
%     SIGnoaam_6597 = ('p < 0.1');
% else
%     SIGnoaam_6597 = ('p > 0.1');
% end
% strnoaam_6597 = mat2str(trendnoaam_6597); strErnoaam_6597 = mat2str(erronoaam_6597); 
% disp(['NOAA 6597 Tend: ',strnoaam_6597,' \pm ',strErnoaam_6597,'(ms^-^1/decade)']); disp(SIGnoaam_6597)
% clear strnoaam_6597; clear strErnoaam_6597; clear SIGnoaam_6597; clear p_valuenoaam_6597; 
% clear trendnoaam_6597; clear erronoaam_6597; clear lsqnoaam_6597; clear yonoaamy_6597;
% 
% lsqeram_6597 = fitlm(meses_6597,sam_era_6597);
% yoeram_6597 = lsqeram_6597.Coefficients.Estimate(1);
% trenderam_6597 = lsqeram_6597.Coefficients.Estimate(2);
% erroeram_6597 = lsqeram_6597.Coefficients.SE(2);
% lineeram_6597 = yoeram_6597+trenderam_6597*meses_6597; 
% trenderam_6597 = (trenderam_6597*120); erroeram_6597 = erroeram_6597*120;
% p_valueeram_6597 = lsqeram_6597.Coefficients.pValue(2); 
% 
% if p_valueeram_6597 < 0.01
%     SIGeram_6597 = ('p < 0.01');
% elseif p_valueeram_6597 < 0.05
%     SIGeram_6597 = ('p < 0.05');
% elseif p_valueeram_6597 < 0.1
%     SIGeram_6597 = ('p < 0.1');
% else
%     SIGeram_6597 = ('p > 0.1');
% end
% streram_6597 = mat2str(trenderam_6597); strEreram_6597 = mat2str(erroeram_6597); 
% disp(['ERA 6597 Tend: ',streram_6597,' \pm ',strEreram_6597]); disp(SIGeram_6597)
% clear streram_6597; clear strEreram_6597; clear SIGeram_6597; clear p_valueeram_6597; 
% clear trenderam_6597; clear erroeram_6597; clear lsqeram_6597; clear yoeram_6597;
% 


%%
% figure('color',[1 1 1],'position',[108 305 900 500]) 
% subplot(2,1,1); plot(tempo,sam_nm,'color', [0.87 0.87 0.87],'linewidth',1);
% hold on; plot(tanos,sam_n,'k','linewidth',1.3);
% plot(tempo,linenoaam,'r','linewidth',1.7) 
% title('SAM (NOAA20CR)'); xlabel('Time');
% datetick('x',10,'keepticks'); xlim([tanos(1) tanos(end)]); ylim([-5 5])
% plot([tempo(1) tempo(end)],[0 0],'Color',[0.6 0.6 0.6],'LineStyle','--')
% text(tanos(2),4,'Trend: 5.02e-04 \pm 9.40-05 (ms^-^1/decade); p < 0.01')
% plot(tempo_6597,linenoaam_6597,'c','linewidth',1.5)
% plot(tempo_8298,linenoaam_8298,'linewidth',1.5)
% plot(tempo_5800,linenoaam_5800,'linewidth',1.5)
% 
% subplot(2,1,2); plot(tempo,sam_em,'color',[0.87 0.87 0.87], 'linewidth',1);
% hold on; plot(tanos,sam_e,'k','linewidth',1.3);
% plot(tanos,lineera,'r','linewidth',1.7)
% title('SAM (ERA20C)'); xlabel('Time');
% datetick('x',10,'keepticks'); xlim([tanos(1) tanos(end)]); ylim([-5 5])
% plot([tempo(1) tempo(end)],[0 0],'Color',[0.6 0.6 0.6],'LineStyle','--')
% text(tanos(2),4,'Trend: 1.23e-03 \pm 8.69e-05 (ms^-^1/decade); p < 0.01')
% plot(tempo_8298,lineeram_8298,'linewidth',1.5)
% plot(tempo_6597,lineeram_6597,'c','linewidth',1.5)
% plot(tempo_5800,lineeram_5800,'linewidth',1.5)