tic
clear;
load('G:\China C neutrality\Data\ID_Pro_CN.mat'); % ID_Pro_CN(4800x1950)
EZZ_CN=imread('G:\China C neutrality\Data\EZZ_CN_all.tif'); % 5-55°N，73-138°E
EZZ_CN = double(EZZ_CN);
% Special Economic zones:  
% 48954 is Japan/China/Taiwan overlap area;  
% 8321 is Taiwan's special economic zone;  
% 8486 is China's special economic zone;  
% 49003 is the overlapping area of the South China Sea;  
EZZ_CN(EZZ_CN~=8486)=0;

SRTM30_CN_sea = zeros(6000,7800);
SRTM30_CN=imread('G:\China C neutrality\Data\SRTM30_CN120.tif');
SRTM30_CN = double(SRTM30_CN);
[m,n]=find(EZZ_CN==8486 & SRTM30_CN<=0);
SRTM30_CN_sea(sub2ind(size(SRTM30_CN_sea), m, n))= SRTM30_CN(sub2ind(size(SRTM30_CN), m, n));

WDPA_CN_sea = zeros(6000,7800);
WDPA_CN = imread('G:\China C neutrality\Data\WDPA_CN120.tif');
WDPA_CN=double(WDPA_CN);
[m,n]=find(WDPA_CN~=0 & EZZ_CN==8486 );
WDPA_CN_sea(sub2ind(size(WDPA_CN_sea), m, n))= WDPA_CN(sub2ind(size(WDPA_CN), m, n));

nature_reserve_sea = zeros(6000,7800);
nature_reserve = imread('G:\China C neutrality\Data\Nature_reserve120.tif'); 
nature_reserve=double(nature_reserve);
[m,n]=find(nature_reserve==47 | nature_reserve==48 | nature_reserve==49 | nature_reserve==50); % Marine ecological function reserve
nature_reserve_sea(sub2ind(size(nature_reserve_sea), m, n))= nature_reserve(sub2ind(size(nature_reserve), m, n));

Marine_reserve_sea = zeros(6000,7800);
[m,n]=find(nature_reserve_sea~=0 | WDPA_CN_sea~=0);
Marine_reserve_sea(sub2ind(size(Marine_reserve_sea), m, n))= 1;

load('G:\China C neutrality\Data\SO2_mean.mat');  % SO2 anthropogenic emissions  kg/m2/s; Spatial:0.5 °x 0.625 °
for i = 1:100
    for j = 1:104
        SO2_mean120(i*60-59:i*60,j*75-74:j*75)=SO2_mean(i,j);
    end
    i
end


S_offshorewind1 = zeros(6000,7800);
[m,n]=find(EZZ_CN==8486 & SRTM30_CN_sea>=-60 & SRTM30_CN_sea<=0 & Marine_reserve_sea~=1);
S_offshorewind1(sub2ind(size(S_offshorewind1), m, n))= 1;
[m,n]=find(EZZ_CN==8486 & SRTM30_CN_sea>=-60 & SRTM30_CN_sea<=0 & Marine_reserve_sea~=1 & SO2_mean120>10^(-11));
S_offshorewind1(sub2ind(size(S_offshorewind1), m, n))= 0.8;

S_offshorewind = S_offshorewind1(1:4800,:);
save('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\S_offshorewind.dat','S_offshorewind');

S_offshorewind2 = S_offshorewind(1201:4800,32*120+1:57*120); % 15°-45°N， 105°-13°E
save('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\S_offshorewind2.dat','S_offshorewind2');
