tic
clear
discount=0.05; % per year
lifetime_power = 25; % year
lifetime_inverter=10; % renewed per 10 years
OMratio_PV = 0.01;
OMratio_wind=0.03;
discountinverter=0;
for t=1:floor(lifetime_power/lifetime_inverter)
    discountinverter=discountinverter+1/(1+discount)^((t-1)*lifetime_inverter);
end

discount40yr=zeros(4,1);
for t=1:lifetime_power
    discount40yr(4,1)=discount40yr(4,1)+1/(1+discount)^(t-1);
end

disrate_PV = ones(4,1)+OMratio_PV.*discount40yr;
disrate_inverter_PV = discountinverter.*ones(4,1)+OMratio_PV.*discount40yr;
disrate_wind = ones(4,1)+OMratio_wind.*discount40yr;

degration_PV = 0.00;
degration_onshorewind = 0.00;
degration_offshorewind = 0.0; 
degrat40yr_PV=zeros(4,1);
degrat40yr_onshorewind=zeros(4,1);
degrat40yr_offshorewind=zeros(4,1);

for t=1:lifetime_power
    degrat40yr_PV(4,1)=degrat40yr_PV(4,1)+(1-degration_PV)^(t-1)/(1+discount)^(t-1);
    degrat40yr_onshorewind(4,1)=degrat40yr_onshorewind(4,1)+(1-degration_onshorewind)^(t-1)/(1+discount)^(t-1);
    degrat40yr_offshorewind(4,1)=degrat40yr_offshorewind(4,1)+(1-degration_offshorewind)^(t-1)/(1+discount)^(t-1);
end

load('H:\China C neutrality\PV_power potential\ANS_PV1\optpowerunit_PV.mat'); % 
load('H:\China C neutrality\PV_power potential\ANS_PV1\powerunit_IX_PV.mat'); % 
load('H:\China C neutrality\PV_power potential\ANS_PV1\tranmission_lines_IX.mat');  % lines_IX
lines_IX(size(optpowerunit_PV,1)+1:end,:)=[];
lines_IX_PV = lines_IX;
load('H:\China C neutrality\PV_power potential\ANS_PV1\powerunit_num_IX_PV.mat');  % 
load('H:\China C neutrality\PV_power potential\ANS_PV1\unitid_lcoe.dat','-mat'); 
unitid_lcoe_PV = unitid_lcoe;
optpowerunit_PV(:,35) = 1;
optpowerunit_PV(:,40) = powerunit_IX_PV; 

load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\optpowerunit_onshorewind.mat'); % 
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit_IX_onshorewind.mat'); % 
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\tranmission_lines_IX.mat');  % lines_IX
lines_IX(size(optpowerunit_onshorewind,1)+1:end,:)=[];
lines_IX_onshorewind = lines_IX;
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit_num_IX_onshorewind.mat');  % 
optpowerunit_onshorewind(:,35) = 2;
optpowerunit_onshorewind(:,40) = powerunit_IX_onshorewind;

load('H:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\optpowerunit_offshorewind.mat'); % 
load('H:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\powerunit_IX_offshorewind.mat'); % 
load('H:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\tranmission_lines_IX.mat');  % lines_IX
lines_IX(size(powerunit_IX_offshorewind,1)+1:end,:)=[];
lines_IX_offshorewind = lines_IX;

optpowerunit_offshorewind(:,20)=optpowerunit_offshorewind(:,8);
optpowerunit_offshorewind(:,30)=optpowerunit_offshorewind(:,3)/1000; %MW
optpowerunit_offshorewind(:,35) = 3;
optpowerunit_offshorewind(:,40) = powerunit_IX_offshorewind;

optpowerunit = [optpowerunit_PV;optpowerunit_onshorewind;optpowerunit_offshorewind];
lines_IX = [lines_IX_PV(:,1:15);lines_IX_onshorewind;lines_IX_offshorewind];
[B,IX]=sort(optpowerunit(:,20),1);
numpowerunit = size(optpowerunit,1);
for i=1:numpowerunit
    i2=IX(i);
    powerunit_IX(i,1)=i2;
    optpowerunit_IX(i,1:40)=optpowerunit(i2,1:40); % lat lon
    lines_IX_IX(i,1:15)=lines_IX(i2,1:15); % lat lon
end

%%
load('H:\China C neutrality\ANS\CO2_mecha.mat'); % USD2019/kWh
load('H:\China C neutrality\ANS\B_mecha.mat'); % USD2019/kWh
load('H:\China C neutrality\ANS\LCOEE_mecha.mat'); % USD2019/kWh
CO2_mecha = CO2all_c_utilize_trans_storage;
B_mecha = B_utilize_trans_storage;
LCOE_mecha = LCOEE_all_utilize_trans_storage;

load('H:\China C neutrality\ANS\CO2_Battery.mat'); % USD2019/kWh
load('H:\China C neutrality\ANS\B_Battery.mat'); % USD2019/kWh
load('H:\China C neutrality\ANS\LCOEE_Battery.mat'); % USD2019/kWh
CO2_Battery = CO2all_c_utilize_trans_storage;
B_Battery = B_utilize_trans_storage;
LCOE_Battery = LCOEE_all_utilize_trans_storage;

CO2=zeros(3844,1);
B=zeros(3844,1);
for i = 1:3844
    if LCOE_mecha(i,1)<LCOE_Battery(i,1)
        cho(i,1)=1; %机械储能
        B(i,1) = B_mecha(i,1);
        CO2(i,1) = CO2_mecha(i,1);
    end
    if LCOE_mecha(i,1)>LCOE_Battery(i,1)
        cho(i,1)=-1; % 化学储能
        B(i,1) = B_Battery(i,1);
        CO2(i,1) = CO2_Battery(i,1);
    end
    if LCOE_mecha(i,1)==LCOE_Battery(i,1)
        cho(i,1)=0; % 无需储能
        B(i,1) = B_Battery(i,1);
        CO2(i,1) = CO2_Battery(i,1);
    end
end
save('H:\China C neutrality\ANS\cho.mat','cho'); % 

[B_IX,IX]=sort(B);
CO2_IX = CO2(IX,1);
CO2_Battery_IX = CO2_Battery(IX,1);
B_Battery_IX = B_Battery(IX,1);
CO2_mecha_IX = CO2_mecha(IX,1);
B_mecha_IX = B_mecha(IX,1);

%%
CP111=optpowerunit_IX(:,30)/10^6; %TW
CP11122 = CP111(IX,1);
CP_cumsum_all=cumsum(CP11122);
[~,Index] = min(abs(CP_cumsum_all-8));

ID_Plant(:,1) =optpowerunit_IX(:,40);
ID_Plant(:,2) =optpowerunit_IX(:,35); % type ,1-3
ID_Plant_IX = ID_Plant(IX,:);
ID_Plant_IX_8TW = ID_Plant_IX(1:Index,:);
save('H:\China C neutrality\ANS\ID_Plant_IX_8TW.mat','ID_Plant_IX_8TW'); % 
