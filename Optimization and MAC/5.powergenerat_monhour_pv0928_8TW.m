tic
clear;

load('H:\China C neutrality\PV_power potential\ANS_PV1\tranmission_lines.dat','-mat');  % lines
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
lines(numlines+1:end,:)=[];

load('H:\China C neutrality\PV_power potential\ANS_PV1\optpowerunit_PV.mat'); % 
load('H:\China C neutrality\PV_power potential\ANS_PV1\powerunit_IX_PV.mat'); % 
load('H:\China C neutrality\PV_power potential\ANS_PV1\tranmission_lines_IX.mat');  % lines_IX
lines_IX(size(optpowerunit_PV,1)+1:end,:)=[];
load('H:\China C neutrality\PV_power potential\ANS_PV1\powerunit_num_IX_PV.mat');  % 
optpowerunit_PV(:,35) = 1;
optpowerunit_PV(:,40) = powerunit_IX_PV; % power plant type ID
[B,IX]=sort(optpowerunit_PV(:,40));
lines_IX_PV_IX = lines_IX(IX,:);
pro_ix = lines_IX_PV_IX(:,10);

%
load('H:\China C neutrality\Data\powergenerat_monhour_pv0928.mat')  % TWh/h  
load('H:\China C neutrality\ANS\ID_Plant_IX_8TW.mat'); % 
[m,n] = find(ID_Plant_IX_8TW(:,2)==1);
id = ID_Plant_IX_8TW(m,1); % 太阳能电厂的序号
for i = 1:size(id,1)
    id(i,2) = pro_ix(id(i,1),1);
end

for i = 1:34
    [m,n] = find(id(:,2)==i);
    powergenerat_monhour_pv_8TW_pro(:,i) = sum(powergenerat_monhour_pv(:,id(m,1)),2);
end

save('H:\China C neutrality\ANS\powergenerat_monhour_pv_8TW_pro.mat','powergenerat_monhour_pv_8TW_pro')  % 288*34 % TWh/h  
%%
tic
clear;

load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\optpowerunit_onshorewind.mat'); % 
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit_IX_onshorewind.mat'); % 
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\tranmission_lines_IX.mat');  % lines_IX
lines_IX(size(optpowerunit_onshorewind,1)+1:end,:)=[];
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit_num_IX_onshorewind.mat');  % 
optpowerunit_onshorewind(:,35) = 2;
optpowerunit_onshorewind(:,40) = powerunit_IX_onshorewind; % power plant type ID
[B,IX]=sort(optpowerunit_onshorewind(:,40));
lines_IX_onshorewind = lines_IX(IX,:);
pro_ix = lines_IX_onshorewind(:,10);

%
load('H:\China C neutrality\Data\powergenerat_monhour_onshorewind0928.mat')   % TWh/h  
load('H:\China C neutrality\ANS\ID_Plant_IX_8TW.mat'); % 
[m,n] = find(ID_Plant_IX_8TW(:,2)==2);
id = ID_Plant_IX_8TW(m,1); % onshorewind电厂的序号
for i = 1:size(id,1)
    id(i,2) = pro_ix(id(i,1),1);
end

for i = 1:34
    [m,n] = find(id(:,2)==i);
    powergenerat_monhour_onshorewind_8TW_pro(:,i) = sum(powergenerat_monhour_onshorewind(:,id(m,1)),2);
end
save('H:\China C neutrality\ANS\powergenerat_monhour_onshorewind_8TW_pro.mat','powergenerat_monhour_onshorewind_8TW_pro')  % 288*34 % TWh/h  

%%
tic
clear;
load('H:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\optpowerunit_offshorewind.mat'); % 
load('H:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\powerunit_IX_offshorewind.mat'); % 
load('H:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\tranmission_lines_IX.mat');  % lines_IX
lines_IX(size(powerunit_IX_offshorewind,1)+1:end,:)=[];
lines_IX_offshorewind = lines_IX;
load('H:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\powerunit_num_IX_offshorewind.mat');  % 
powerunit_num_IX_offshorewind(:,4) = 1;
optpowerunit_offshorewind(:,20)=optpowerunit_offshorewind(:,8);
optpowerunit_offshorewind(:,30)=optpowerunit_offshorewind(:,3)/1000; %MW
optpowerunit_offshorewind(:,35) = 3;
optpowerunit_offshorewind(:,40) = powerunit_IX_offshorewind; % power plant type ID

[B,IX]=sort(optpowerunit_offshorewind(:,40));
lines_IX_offshorewind = lines_IX(IX,:);
pro_ix = lines_IX_offshorewind(:,10);

load('H:\China C neutrality\Data\powergenerat_monhour_offshorewind0928.mat')  % TWh/h  
load('H:\China C neutrality\ANS\ID_Plant_IX_8TW.mat'); % 
[m,n] = find(ID_Plant_IX_8TW(:,2)==3);
id = ID_Plant_IX_8TW(m,1); % onshorewind电厂的序号
for i = 1:size(id,1)
    id(i,2) = pro_ix(id(i,1),1);
end

for i = 1:34
    [m,n] = find(id(:,2)==i);
    powergenerat_monhour_offshorewind_8TW_pro(:,i) = sum(powergenerat_monhour_offshorewind(:,id(m,1)),2);
end

save('H:\China C neutrality\ANS\powergenerat_monhour_offshorewind_8TW_pro.mat','powergenerat_monhour_offshorewind_8TW_pro')  % 288*34 % TWh/h  
