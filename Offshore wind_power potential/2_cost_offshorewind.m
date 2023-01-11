tic
clear;
load('H:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\S_offshorewind.dat','-mat'); 
load('H:\China C neutrality\Data\SRTM30_CN120.mat');
SRTM30_CN = double(SRTM30_CN);
SRTM30_CN = SRTM30_CN(1:4800,:);
SRTM30_CN1 = zeros(4800,7800);
[m,n]=find(S_offshorewind~=0 & SRTM30_CN<=0);
SRTM30_CN1(sub2ind(size(SRTM30_CN1), m, n))= SRTM30_CN(sub2ind(size(SRTM30_CN), m, n));


load('H:\China C neutrality\Data\dist_S_offshore.mat'); 
% 2: m; 3: n; 
% 4：Lat；5：Lon；
% 6：NEAR_DIST(m)
% 7：Pro_ID
dist_offshore = zeros(40*30,65*30);
pro_offshore = zeros(40*30,65*30);
dist_offshore(sub2ind(size(dist_offshore), dist_S_offshore(:,2), dist_S_offshore(:,3)))= dist_S_offshore(:,6)/1000;
pro_offshore(sub2ind(size(dist_offshore), dist_S_offshore(:,2), dist_S_offshore(:,3)))= dist_S_offshore(:,7);

for i = 1:40*30
    for j = 1:65*30
        dist_offshore120(i*4-3:i*4,j*4-3:j*4)=dist_offshore(i,j); % km
        pro_offshore120(i*4-3:i*4,j*4-3:j*4)=pro_offshore(i,j);
    end
    i
end

C0 = 2000/1000; % $/W
C = C0.*(0.0084.*(-SRTM30_CN1)+0.8368).*(0.0057.*dist_offshore120+0.7714);
cost_offshorewind = zeros(4800,7800);
[m,n]=find(S_offshorewind~=0 & SRTM30_CN<=0);
cost_offshorewind(sub2ind(size(cost_offshorewind), m, n))= C(sub2ind(size(C), m, n));

save('H:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\cost_offshorewind.mat','cost_offshorewind'); % $/W
save('H:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\dist_offshore120.mat','dist_offshore120'); % km
save('H:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\pro_offshore120.mat','pro_offshore120'); % pro ID
