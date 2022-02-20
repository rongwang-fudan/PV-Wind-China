tic
clear;
load('G:\China C neutrality\Data\Wind_speed_100m_mean120.mat'); % m/s
Wind_speed_100m_mean120(find(isnan(Wind_speed_100m_mean120)==1))=0;

load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\dist_offshore120.mat'); % km
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\pro_offshore120.mat'); % pro ID
% 3:Fujian; 5:Guangdong; 6: Guangxi; 8:Hainan; 9:Hebei; 17:Jiangsu; 
% 18:Liaoning; 25: Shandong; 26:Shanghai; 28: Tianjin; 32: Yunnan; 33:Zhejiang; 
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\cost_offshorewind.mat'); % $/W
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\S_offshorewind.dat','-mat'); 
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\ID_dist_offshore120.mat'); % $/W
SRTM30_CN=imread('G:\China C neutrality\Data\SRTM30_CN120.tif');
SRTM30_CN = double(SRTM30_CN);
SRTM30_CN = SRTM30_CN(1:4800,:);
SRTM30_CN120 = zeros(4800,7800);
[m,n]=find(S_offshorewind~=0 & SRTM30_CN<=0);
SRTM30_CN120(sub2ind(size(SRTM30_CN120), m, n))= SRTM30_CN(sub2ind(size(SRTM30_CN), m, n));

D = 164; %m
Pwp = 8*10^6;%W
CP_unit = Pwp/(7*D*7*D); %W/m2
load('G:\China C neutrality\Data\gridarea120.mat'); 
grid_area = gridarea120 *ones(1,65*120)*10^6; % m2
S_offshorewind=grid_area.*S_offshorewind;% m2
CP = CP_unit.*S_offshorewind/1000; % kW
Cost = cost_offshorewind.*CP*1000/10^6; % million dollar
% Ele_coef120
Depth = SRTM30_CN120.*(-1);
Ele_loss =(2.07+(0.073*dist_offshore120)+(-0.0016*dist_offshore120.^2 )+(0.000017*dist_offshore120.^3 )+(-0.000000086*dist_offshore120.^4 )+(0.000000000157*dist_offshore120.^5 )+0.0015*Depth+(-0.0000047*Depth.^2 )+(0.0000000082*Depth.^3 )+(-0.0000000000041*Depth.^4 ))/100;
Ele_coef = ones(4800,7800)-Ele_loss;
Ele_coef120 = zeros(4800,7800);
[m,n]=find(S_offshorewind~=0 & SRTM30_CN<=0);
Ele_coef120(sub2ind(size(Ele_coef120), m, n))= Ele_coef(sub2ind(size(Ele_coef), m, n));

UTI_coef=0.95;
ARR_coef=0.90;
Other_coef=0.98;

load('G:\China C neutrality\Data\CF_offshore.mat'); % 
Ph = CP.*CF_offshore.*UTI_coef.*ARR_coef.*Ele_coef120.*Other_coef.*8760; % kwh/year
save('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\offshorewind_CP_all.mat','CP'); % kW
save('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\offshorewind_power_all.mat','Ph'); % kwh/year


tt=1;
for i=[3 5 6 8 9 17 18 25 26 28 32 33]
    ID_dist_offshore120_pro = zeros(4800,7800);
    [m,n]=find(S_offshorewind~=0 & SRTM30_CN<=0 & pro_offshore120==i);
    ID_dist_offshore120_pro(sub2ind(size(ID_dist_offshore120_pro), m, n))= ID_dist_offshore120(sub2ind(size(ID_dist_offshore120), m, n));
    for j=1:max(max(ID_dist_offshore120_pro))
        [m,n]=find(ID_dist_offshore120_pro==j);
        costunits_offshorewind(tt,1)= i; % pro
        costunits_offshorewind(tt,2)= j; % Powerunit_ID
        costunits_offshorewind(tt,3)= sum(sum(CP(sub2ind(size(CP), m, n)))); % CP,kW
        costunits_offshorewind(tt,4)= sum(sum(Cost(sub2ind(size(Cost), m, n))));% Cost,million dollar
        costunits_offshorewind(tt,5)= sum(sum(Ph(sub2ind(size(Ph), m, n)))); % kwh/year
        costunits_offshorewind(tt,6)= sum(sum(S_offshorewind(sub2ind(size(S_offshorewind), m, n))))/10^6; % km2
        costunits_offshorewind(tt,7)= sum(sum(Wind_speed_100m_mean120(sub2ind(size(Wind_speed_100m_mean120), m, n))))/size(m,1); % mean wind speed, m/s
        tt=tt+1;
    end
    i
end
save('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\costunits_offshorewind.dat','costunits_offshorewind');
    