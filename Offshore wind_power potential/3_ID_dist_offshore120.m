tic
clear;
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\dist_offshore120.mat'); % km
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\pro_offshore120.mat'); % pro ID
% 3:Fujian; 5:Guangdong; 6: Guangxi; 8:Hainan; 9:Hebei; 17:Jiangsu; 
% 18:Liaoning; 25: Shandong; 26:Shanghai; 28: Tianjin; 32: Yunnan; 33:Zhejiang; 
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\cost_offshorewind.mat'); % $/W

dist_offshore120(dist_offshore120==0)=1;
dist_offshore120=ceil(dist_offshore120);

load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\S_offshorewind.dat','-mat'); 
SRTM30_CN=imread('G:\China C neutrality\Data\SRTM30_CN120.tif');
SRTM30_CN = double(SRTM30_CN);
SRTM30_CN = SRTM30_CN(1:4800,:);
ID_dist_offshore120 = zeros(4800,7800);
[m,n]=find(S_offshorewind~=0 & SRTM30_CN<=0);
ID_dist_offshore120(sub2ind(size(ID_dist_offshore120), m, n))= dist_offshore120(sub2ind(size(dist_offshore120), m, n));
save('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\ID_dist_offshore120.mat','ID_dist_offshore120'); % power ID
