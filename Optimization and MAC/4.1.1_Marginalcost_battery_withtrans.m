tic
clear;
Cost_mechanical = zeros(3844,1); %

load('H:\China C neutrality\PV_power potential\ANS_PV1\tranmission_lines.dat','-mat');  % lines
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
lines(numlines+1:end,:)=[];
load('H:\China C neutrality\ANS\unitmin_module.mat'); % The number of PV,onshore wind and offshore wind power plants at each stage of construction at minimum cost
load('H:\China C neutrality\ANS\utilize_ratio2060_trans_plant_alone.mat')  % power use efficiency with inter-regional electric transport but without storage
load('H:\China C neutrality\ANS\utilize_ratio2060_trans_storage_plant_alone.mat')  % power use efficiency with inter-regional electric transport and storage
load('H:\China C neutrality\ANS\storage_max_plant.mat')   % Maximum charging power, TWh/h  
load('H:\China C neutrality\ANS\storage_year_plant.mat')  % Average annual electricity storage, TWh/year
load('H:\China C neutrality\ANS\cost_trans_IX2.mat')% transmission cost
load('H:\China C neutrality\ANS\CP_trans_IX2.mat')% 

lifetime_s = 15; % year
lifetime_power = 25; % year
Ip = storage_max_plant/(0.99*0.85*0.983*0.967)*10^9; % kW  power capacity
Ie = storage_year_plant/(0.99*0.85*0.983*0.967)*10^9*lifetime_s/6000; % kWh  energy capacity
for i = 1:3844
    if Ip(i,1)>Ie(i,1)
        Ie(i,1) = Ip(i,1);
    end
end

p_dis = storage_year_plant*10^9/(0.983*0.967); % annual discharge, kWh/year 
p_char = storage_year_plant*10^9/(0.99*0.85*0.983*0.967); % annual charge, kWh/year 
Co = 1.5/1000; % $/kWh
Cp=595.73; % $/kW
Ce = 345; % 2018$/kWh
discount=0.05; % per year
discount1yr=0;
for t=1:lifetime_power
    discount1yr=discount1yr+1/(1+discount)^(t-1);
end
discount2yr=0;
for t=16
    discount2yr=discount2yr+1/(1+discount)^(t-1);
end
r_y = (lifetime_power-lifetime_s)/lifetime_s;

unitmin = [0;unitmin];
rIp=[595.73, 374.45, 327.22,280.78]./595.73;
rIe=[345,198, 174, 149]./345;
for i = 1:4
    learnprice_sto_stage_p(unitmin(i)+1:unitmin(i+1),:) = rIp(i);
    learnprice_sto_stage_e(unitmin(i)+1:unitmin(i+1),:) = rIe(i);
end

Cost_storage1 = (Ce.*Ie.*learnprice_sto_stage_e+Cp.*Ip.*learnprice_sto_stage_p)/10^6+(Ce.*Ie.*learnprice_sto_stage_e+Cp.*Ip.*learnprice_sto_stage_p)/10^6*discount2yr*r_y+Co*(p_dis+p_char)/10^6*discount1yr; %million $
clear discount1yr
Cost_storage=Cost_storage1;
clear Cost_storage1

%% trans
LR_PV1 = 0; %0.18; % 0.37; %  0.37 PV module
LR_PV2 = 0.18; %0.18; % 0.37; %  0.37 PV module
CP0_PV = sum(sum(CP_trans_IX))/10; % MW
CP_PV = zeros(size(Ie,1),1);
CP_PV(1,1) = CP0_PV; %kw
for i = 1: size(CP_trans_IX,2)
        CP_PV(i+1,1) =  CP_PV(i,1) + CP_trans_IX(i);
end
r_pv = CP_PV./CP0_PV;
for i = 1:size(r_pv,1)
    if CP_PV(i,1)<=CP0_PV
        LR_trans(i,1) = LR_PV1;
        learnprice_trans(i,1) = (CP_PV(i,1)./CP0_PV).^(log(1- LR_trans(i,1))/log(2));
    else
        LR_trans(i,1) = LR_PV2;
        learnprice_trans(i,1) = (CP_PV(i,1)./CP0_PV).^(log(1- LR_trans(i,1))/log(2));      
    end
end

for i = 1:4
    learnprice_trans2(unitmin(i)+1:unitmin(i+1),:) = learnprice_trans(unitmin(i)+1,1);
end
learnprice_trans = learnprice_trans2;
clear learnprice_trans2


%%
CO2_C=0.2727;
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
CF=sum(optpowerunit_IX(:,1))./sum(optpowerunit_IX(:,30))/8760*10^6;
[B,IX]=sort(optpowerunit_IX(:,20),1);

%% caculate power use efficiency
utilize_ra = (utilize_ratio2060_trans_plant_alone.*optpowerunit_IX(:,1)+storage_year_plant/(0.99*0.85*0.983*0.967)*(0.99*0.85*0.983*0.967))./optpowerunit_IX(:,1);
utilize_ratio2060_trans_plant_alone = utilize_ra;
utilize_ratio2060_trans_storage_plant_alone = utilize_ra;
utilize_ratio2060_2 = utilize_ra;
clear utilize_ra
%%
line_IX_all = [lines(:,1:15);lines_IX_IX];
[m,n] = find(line_IX_all(:,9)==0);
line_IX_all(m,9) = 1362; % offshorewind, Shenzhen

% PV
LR_PV1 = 0.2; 
LR_PV2 = 0.18;
CP0_PV = 253000;  % cumulative capacity potential in China in 2020, MW
CP_PV = zeros(numpowerunit,1);
CP_PV(1,1) = CP0_PV;
for i = 1: numpowerunit
    if optpowerunit_IX(i,35)==1
        CP_PV(i+1,1) =  CP_PV(i,1) + optpowerunit_IX(i,30);
    else
        CP_PV(i+1,1) =  CP_PV(i,1);
    end    
end
CP0_PV_cum = CP_PV(2:end,1)-CP0_PV;
r_pv = CP_PV./CP0_PV;
for i = 1:size(r_pv,1)
        LR_PV(i,1) = LR_PV1;
        learnprice_PV(i,1) = (CP_PV(i,1)./CP0_PV).^(log(1- LR_PV(i,1))/log(2));
        LR_PV_2(i,1) = LR_PV2;
        learnprice_PV_2(i,1) = (CP_PV(i,1)./CP0_PV).^(log(1- LR_PV_2(i,1))/log(2));        
end

for i = 1:size(CP_PV,1)-1
    CP_PV22(i,1) = CP_PV(i+1,1)-CP_PV(1,1);
    CP_PV22(i,2) = (CP_PV(i+1,1)-CP_PV(1,1))/(CP_PV(end,1)-CP_PV(1,1));
end
cost_PV = zeros(numpowerunit,19);
cost_PV_alone = zeros(numpowerunit,30);
cost_PV_trans = zeros(numpowerunit,1);
CO2_PV_alone = zeros(numpowerunit,3);
if optpowerunit_IX(1,35)==1
    cost_PV(1,11:19) =  optpowerunit_IX(1,11:19);
    cost_PV_alone(1,1:30) = optpowerunit_IX(1,1:30);
    cost_PV_trans(1,1) = cost_trans_IX(1,1);
    CO2_PV_alone(1,1:3)  = optpowerunit_IX(1,8:10);
end
for i = 2: numpowerunit
    if optpowerunit_IX(i,35)==1
        cost_PV(i,11:19) =  cost_PV(i-1,11:19) + optpowerunit_IX(i,11:19);
        cost_PV_alone(i,1:30) =  optpowerunit_IX(i,1:30);
        cost_PV_trans(i,1) = cost_trans_IX(i,1);
        CO2_PV_alone(i,1:3)  = optpowerunit_IX(i,8:10);
    else
        cost_PV(i,11:19) =  cost_PV(i-1,11:19);
        cost_PV_alone(i,1:30) =  0;
        cost_PV_trans(i,1) = 0;
        CO2_PV_alone(i,1:3)  = 0;
    end    
end

% onshorewind
LR_onshorewind1= 0.074; 
LR_onshorewind2= 0.18;
CP0_onshorewind = 272010;% MW
CP_onshorewind = zeros(numpowerunit,1);
CP_onshorewind(1,1) = CP0_onshorewind;
for i = 1: numpowerunit
    if optpowerunit_IX(i,35)==2
        CP_onshorewind(i+1,1) =  CP_onshorewind(i,1) + optpowerunit_IX(i,30);
    else
        CP_onshorewind(i+1,1) =  CP_onshorewind(i,1);
    end    
end
CP0_onshorewind_cum = CP_onshorewind(2:end,1)-CP0_onshorewind;
r_onshorewind = CP_onshorewind./CP0_onshorewind;
for i = 1:size(r_onshorewind,1)
        LR_onshorewind(i,1) = LR_onshorewind1;
        learnprice_onshorewind(i,1) = (CP_onshorewind(i,1)./CP0_onshorewind).^(log(1- LR_onshorewind(i,1))/log(2));
        LR_onshorewind_2(i,1) = LR_onshorewind2;
        learnprice_onshorewind_2(i,1) = (CP_onshorewind(i,1)./CP0_onshorewind).^(log(1- LR_onshorewind_2(i,1))/log(2));
end

for i = 1:size(CP_onshorewind,1)-1
    CP_onshorewind22(i,1) = CP_onshorewind(i+1,1)-CP_onshorewind(1,1);
    CP_onshorewind22(i,2) = (CP_onshorewind(i+1,1)-CP_onshorewind(1,1))/(CP_onshorewind(end,1)-CP_onshorewind(1,1));
end
cost_onshorewind = zeros(numpowerunit,19);
cost_onshorewind_alone = zeros(numpowerunit,30);
cost_onshorewind_trans = zeros(numpowerunit,1);
CO2_onshorewind_alone = zeros(numpowerunit,3);
if optpowerunit_IX(1,35)==2
    cost_onshorewind(1,11:19) =  optpowerunit_IX(1,11:19);
    cost_onshorewind_alone(1,1:30) =  optpowerunit_IX(1,1:30);
    cost_onshorewind_trans(1,1) = cost_trans_IX(1,1);
    CO2_onshorewind_alone(1,1:3) =  optpowerunit_IX(1,8:10);
end
for i = 2: numpowerunit
    if optpowerunit_IX(i,35)==2
        cost_onshorewind(i,11:19) =  cost_onshorewind(i-1,11:19) + optpowerunit_IX(i,11:19);
        cost_onshorewind_alone(i,1:30) =  optpowerunit_IX(i,1:30);
        cost_onshorewind_trans(i,1) = cost_trans_IX(i,1);
        CO2_onshorewind_alone(i,1:3) =  optpowerunit_IX(i,8:10);
    else
        cost_onshorewind(i,11:19) =  cost_onshorewind(i-1,11:19);
        cost_onshorewind_alone(i,1:30) =  0;
        cost_onshorewind_trans(i,1) = 0;
        CO2_onshorewind_alone(i,1:3) =  0;
    end    
end

% offhorewind
LR_offshorewind1 = 0.074;
LR_offshorewind2 = 0.18;
CP0_offshorewind = 8990;
CP_offshorewind = zeros(numpowerunit,1);
CP_offshorewind(1,1) = CP0_offshorewind;
for i = 1: numpowerunit
    if optpowerunit_IX(i,35)==3
        CP_offshorewind(i+1,1) =  CP_offshorewind(i,1) + optpowerunit_IX(i,30);
    else
        CP_offshorewind(i+1,1) =  CP_offshorewind(i,1);
    end    
end
CP0_offshorewind_cum = CP_offshorewind(2:end,1)-CP0_offshorewind;
for i = 1:size(CP_offshorewind,1)-1 
    CP_offshorewind22(i,1) = CP_offshorewind(i+1,1)-CP_offshorewind(1,1);
    CP_offshorewind22(i,2) = (CP_offshorewind(i+1,1)-CP_offshorewind(1,1))/(CP_offshorewind(end,1)-CP_offshorewind(1,1));
end
r_offshorewind = CP_offshorewind./CP0_offshorewind;
CP_offshorewind=CP_offshorewind;

for i = 1:size(r_offshorewind,1)
        LR_offshorewind(i,1) = LR_offshorewind1;
        learnprice_offshorewind(i,1) = (CP_offshorewind(i,1)./CP0_offshorewind).^(log(1- LR_offshorewind(i,1))/log(2));
        LR_offshorewind_2(i,1) = LR_offshorewind2;
        learnprice_offshorewind_2(i,1) = (CP_offshorewind(i,1)./CP0_offshorewind).^(log(1- LR_offshorewind_2(i,1))/log(2));
end

cost_offshorewind = zeros(numpowerunit,10);
cost_offshorewind_alone = zeros(numpowerunit,30);
cost_offshorewind_trans = zeros(numpowerunit,1);
if optpowerunit_IX(1,35)==3
    cost_offshorewind(1,6:10) =  optpowerunit_IX(1,6:10);
    cost_offshorewind_alone(1,1:30) =  optpowerunit_IX(1,1:30);
    cost_offshorewind_trans(1,1) = cost_trans_IX(1,1);
end
for i = 2: numpowerunit
    if optpowerunit_IX(i,35)==3
        cost_offshorewind(i,6:10) =  cost_offshorewind(i-1,6:10) + optpowerunit_IX(i,6:10);
        cost_offshorewind_alone(i,1:30) =  optpowerunit_IX(i,1:30);        
        cost_offshorewind_trans(i,1) = cost_trans_IX(i,1);
    else
        cost_offshorewind(i,6:10) =  cost_offshorewind(i-1,6:10);
        cost_offshorewind_alone(i,1:30) =  0;
        cost_offshorewind_trans(i,1) = 0;
    end    
end

optpowerunit1= zeros(numpowerunit,30);
optpowerunit1(1,30)=optpowerunit_IX(1,30);
for i=2:numpowerunit
    optpowerunit1(i,30)=optpowerunit1(i-1,30)+optpowerunit_IX(i,30); % cumulative capacity
end
optpowerunit_IX(:,36)=optpowerunit1(:,30)/optpowerunit1(end,30);

%%
lcoee_PV = zeros(numpowerunit,1);
lcoee_onshorewind = zeros(numpowerunit,1);
lcoee_offshorewind = zeros(numpowerunit,1);

CP_PV_alone(:,1) = cost_PV_alone(:,30); % Capacity potential, MW
CP_onshorewind_alone(:,1) = cost_onshorewind_alone(:,30); % Capacity potential, MW
CP_offshorewind_alone(:,1) = cost_offshorewind_alone(:,30); % Capacity potential, MW

CPPP=CP_PV_alone;
[m,n]=find(CP_PV_alone~=0);
CPPP(sub2ind(size(CPPP), m, n))= CP_PV_alone(sub2ind(size(CP_PV_alone), m, n));
[m,n]=find(CP_onshorewind_alone~=0);
CPPP(sub2ind(size(CPPP), m, n))= CP_onshorewind_alone(sub2ind(size(CP_onshorewind_alone), m, n));
[m,n]=find(CP_offshorewind_alone~=0);
CPPP(sub2ind(size(CPPP), m, n))= CP_offshorewind_alone(sub2ind(size(CP_offshorewind_alone), m, n));
CP111=CPPP/10^6; %TW

CO2_PV_year_abated(:,1)= cost_PV_alone(:,8)./CO2_C; % Mton CO2
CO2_onshorewind_year_abated(:,1)=cost_onshorewind_alone(:,8)./CO2_C; % Mton CO2
CO2_offshorewind_year_abated(:,1)=cost_offshorewind_alone(:,5)./CO2_C; % Mton CO2
CO2_year_abated = CO2_PV_year_abated+CO2_onshorewind_year_abated+CO2_offshorewind_year_abated;

CO2_PV_year_ls(:,1)= cost_PV_alone(:,9)./CO2_C; % Mton CO2
CO2_onshorewind_year_ls(:,1)= cost_onshorewind_alone(:,9)./CO2_C; % Mton CO2
CO2_offshorewind_year_ls(:,1)= 0; % Mton CO2

% 4  land cost
lcoee_PV(:,4) = cost_PV_alone(:,16)./disrate_PV(4,1) ; % cost million USD 
% 6  abated annual CO2 emission
lcoee_PV(:,6) = cost_PV_alone(:,17); % cost million USD 

% 7  land carbon sink
lcoee_PV(:,7) = cost_PV_alone(:,18); % cost million USD 
% 8  land use change
lcoee_PV(:,8) = cost_PV_alone(:,19); % cost million USD 
% 4  land cost
lcoee_onshorewind(:,4) = cost_onshorewind_alone(:,16)./disrate_wind(4,1) ; % cost million USD 
% 6  abated annual CO2 emission
lcoee_onshorewind(:,6) = cost_onshorewind_alone(:,17); % cost million USD 

% 7  land carbon sink
lcoee_onshorewind(:,7) = cost_onshorewind_alone(:,18); % cost million USD 
% 8  land use change
lcoee_onshorewind(:,8) = cost_onshorewind_alone(:,19); % cost million USD 
% 6  abated annual CO2 emission
lcoee_offshorewind(:,6) = cost_offshorewind_alone(:,7); % cost million USD 

% power generation
phh_PV(:,1) = cost_PV_alone(:,1).*degrat40yr_PV(4,1); 

phh_PV_alone(:,1) = cost_PV_alone(:,1); 
sto_PV_sto2(:,1) = cost_PV_alone(:,1).*(1-utilize_ratio2060_trans_storage_plant_alone); 

% power generation
phh_offshorewind(:,1) = cost_offshorewind_alone(:,1).*degrat40yr_offshorewind(4,1); 

% power generation
phh_onshorewind_utilize_trans_storage(:,1) = cost_onshorewind_alone(:,1).*degrat40yr_onshorewind(4,1).*utilize_ratio2060_trans_storage_plant_alone; 

module_price_offshore = 1.008; % USD/W

%%
for i = 1:4
    learnprice_PV_stage(unitmin(i)+1:unitmin(i+1),:) = learnprice_PV(unitmin(i)+1,1);
    learnprice_PV_stage_2(unitmin(i)+1:unitmin(i+1),:) = learnprice_PV_2(unitmin(i)+1,1);
    learnprice_onshorewind_stage(unitmin(i)+1:unitmin(i+1),:) = learnprice_onshorewind(unitmin(i)+1,1);
    learnprice_onshorewind_stage_2(unitmin(i)+1:unitmin(i+1),:) = learnprice_onshorewind_2(unitmin(i)+1,1);
    learnprice_offshorewind_stage(unitmin(i)+1:unitmin(i+1),:) = learnprice_offshorewind(unitmin(i)+1,1);
    learnprice_offshorewind_stage_2(unitmin(i)+1:unitmin(i+1),:) = learnprice_offshorewind_2(unitmin(i)+1,1);
end
learnprice_PV = learnprice_PV_stage;
learnprice_PV_2 = learnprice_PV_stage_2;
learnprice_onshorewind = learnprice_onshorewind_stage;
learnprice_onshorewind_2 = learnprice_onshorewind_stage_2;
learnprice_offshorewind = learnprice_offshorewind_stage;
learnprice_offshorewind_2 = learnprice_offshorewind_stage_2;
clear learnprice_PV_stage
clear learnprice_PV_stage_2
clear learnprice_onshorewind_stage
clear learnprice_onshorewind_stage_2
clear learnprice_offshorewind_stage
clear learnprice_offshorewind_stage_2


% PV
% 1  fixed system cost
lcoee_PV(:,1) = cost_PV_alone(:,14)./disrate_PV(4,1) .*learnprice_PV + cost_PV_alone(:,15)/disrate_inverter_PV(4,1).*discountinverter .*learnprice_PV; % cost million USD 
% 2  substation
lcoee_PV(:,2) = cost_PV_alone(:,12)./disrate_PV(4,1).*learnprice_PV_2 ; % cost million USD 
% 3  Power plant cable
lcoee_PV(:,3) = cost_PV_alone(:,13)./disrate_PV(4,1).*learnprice_PV_2+cost_PV_trans(:,1).*learnprice_trans;; % cost million USD 
% 5  connection to national grid
lcoee_PV(:,5) = cost_PV_alone(:,11)./disrate_PV(4,1).*learnprice_PV_2; % cost million USD 
% 9  O&M cost
lcoee_PV(:,9) = lcoee_PV(:,4).*(disrate_PV(4,1)-1)+(sum(lcoee_PV(:,1:3),2)+lcoee_PV(:,5)).*(disrate_PV(4,1)-1); % cost million USD 

% onshorewind
% 1  fixed system cost
lcoee_onshorewind(:,1) = cost_onshorewind_alone(:,14)./disrate_wind(4,1).*learnprice_onshorewind; % cost million USD 
% 2  substation
lcoee_onshorewind(:,2) = cost_onshorewind_alone(:,12)./disrate_wind(4,1).*learnprice_onshorewind_2; % cost million USD 
% 3  Power plant cable
lcoee_onshorewind(:,3) = cost_onshorewind_alone(:,13)./disrate_wind(4,1).*learnprice_onshorewind_2+cost_onshorewind_trans(:,1).*learnprice_trans; % cost million USD 
% 5  connection to national grid
lcoee_onshorewind(:,5) = cost_onshorewind_alone(:,11)./disrate_wind(4,1).*learnprice_onshorewind_2; % cost million USD 
% 9  O&M cost
lcoee_onshorewind(:,9) = lcoee_onshorewind(:,4).*(disrate_wind(4,1)-1)+(sum(lcoee_onshorewind(:,1:3),2)+lcoee_onshorewind(:,5)).*(disrate_wind(4,1)-1); % cost million USD 

% offshorewind
% 1  module cost
lcoee_offshorewind(:,1) = module_price_offshore.*CP_offshorewind_alone(:,1).*learnprice_offshorewind; % cost million USD 
% 2  others
lcoee_offshorewind(:,2) = (cost_offshorewind_alone(:,6)./disrate_wind(4,1)-module_price_offshore.*CP_offshorewind_alone(:,1)).*learnprice_offshorewind_2+cost_offshorewind_trans(:,1).*learnprice_trans; % cost million USD 
% 9  O&M cost
lcoee_offshorewind(:,9) = sum(lcoee_offshorewind(:,1:2),2).*(disrate_wind(4,1)-1); % cost million USD 

%%
lcoee_PV_utilize_trans_storage = cost_PV_alone(:,17).*utilize_ratio2060_trans_storage_plant_alone; % cost million USD 
phh_PV_utilize_trans_storage(:,1) = cost_PV_alone(:,1).*degrat40yr_PV(4,1).*utilize_ratio2060_trans_storage_plant_alone; 
lcoee_onshorewind_utilize_trans_storage = cost_onshorewind_alone(:,17).*utilize_ratio2060_trans_storage_plant_alone; % cost million USD 
phh_onshorewind_utilize_trans_storage(:,1) = cost_onshorewind_alone(:,1).*degrat40yr_onshorewind(4,1).*utilize_ratio2060_trans_storage_plant_alone; 
lcoee_offshorewind_utilize_trans_storage = cost_offshorewind_alone(:,7).*utilize_ratio2060_trans_storage_plant_alone; % cost million USD 
phh_offshorewind_utilize_trans_storage(:,1) = cost_offshorewind_alone(:,1).*degrat40yr_offshorewind(4,1).*utilize_ratio2060_trans_storage_plant_alone; 

phhall_utilize_trans_storage = phh_PV_utilize_trans_storage;
lcoeeall_utilize_trans_storage = lcoee_PV;
lcoeeall_utilize_trans_storage(:,6) = lcoee_PV_utilize_trans_storage;
[m,n]=find(phh_onshorewind_utilize_trans_storage~=0);
phhall_utilize_trans_storage(sub2ind(size(phhall_utilize_trans_storage), m, n))= phh_onshorewind_utilize_trans_storage(sub2ind(size(phh_onshorewind_utilize_trans_storage), m, n));
[m,n]=find(lcoee_onshorewind~=0);
lcoeeall_utilize_trans_storage(sub2ind(size(lcoeeall_utilize_trans_storage), m, n))= lcoee_onshorewind(sub2ind(size(lcoee_onshorewind), m, n));
lcoeeall_utilize_trans_storage(sub2ind(size(lcoeeall_utilize_trans_storage), m, 6*ones(size(m,1),1)))= lcoee_onshorewind_utilize_trans_storage(sub2ind(size(lcoee_onshorewind_utilize_trans_storage), m, ones(size(m,1),1)));
[m,n]=find(phh_offshorewind_utilize_trans_storage~=0);
phhall_utilize_trans_storage(sub2ind(size(phhall_utilize_trans_storage), m, n))= phh_offshorewind_utilize_trans_storage(sub2ind(size(phh_offshorewind_utilize_trans_storage), m, n));
[m,n]=find(lcoee_offshorewind~=0);
lcoeeall_utilize_trans_storage(sub2ind(size(lcoeeall_utilize_trans_storage), m, n))= lcoee_offshorewind(sub2ind(size(lcoee_offshorewind), m, n));
lcoeeall_utilize_trans_storage(sub2ind(size(lcoeeall_utilize_trans_storage), m, 6*ones(size(m,1),1)))= lcoee_offshorewind_utilize_trans_storage(sub2ind(size(lcoee_offshorewind_utilize_trans_storage), m, ones(size(m,1),1)));

LCOEE_all_utilize_trans_storage = (sum(lcoeeall_utilize_trans_storage(:,1:9),2)+Cost_storage)./phhall_utilize_trans_storage(:,1)/1000; % LCoE million USD2019/TWh->USD2019/kWh   ;            
lcoeeall_utilize_trans_storage(:,10) = Cost_storage;
save('H:\China C neutrality\ANS\LCOEE_Battery.mat','LCOEE_all_utilize_trans_storage'); % USD2019/kWh

[B,ixx]=sort(LCOEE_all_utilize_trans_storage);
optpowerunit_IX_ix = optpowerunit_IX(ixx,1); 
ixx_2 = zeros(3844,1);
rp = cumsum(optpowerunit_IX_ix(:,1))./sum(optpowerunit_IX_ix(:,1));
[~,Index1] = min(abs(rp-0.9567));
% In_cog(1,1) = Index;
ixx_2(ixx(1:Index1),1)=1;
[~,Index2] = min(abs(rp-0.9974));
ixx_2(ixx(Index1+1:Index2),1)=2;
[~,Index3] = min(abs(rp-1));
ixx_2(ixx(Index2+1:Index3),1)=3;
clear optpowerunit_IX_ix
clear Cost_Aba
clear rp
clear ixx
clear Index1
clear Index2
clear Index3

%% 
fossilfuel_emissionfactor=0.783; 
[m,n]=find(ixx_2==1);
emissionfactor(m,1) = 0.840284356; % coal, kg CO2/kWh
[m,n]=find(ixx_2==2);
emissionfactor(m,1) = 0.455141557; % gas, kg CO2/kWh
[m,n]=find(ixx_2==3);
emissionfactor(m,1) = 0.723540308; % oil, kg CO2/kWh


CO2_PV_all_utilize_trans_storage(:,1)= emissionfactor./fossilfuel_emissionfactor.*(cost_PV_alone(:,8)./CO2_C.*degrat40yr_PV(4,1).*utilize_ratio2060_2+ cost_PV_alone(:,9)./CO2_C.*degrat40yr_PV(4,1)+ cost_PV_alone(:,10)./CO2_C); % Mton CO2
CO2_onshorewind_all_utilize_trans_storage(:,1)=emissionfactor./fossilfuel_emissionfactor.*(cost_onshorewind_alone(:,8)./CO2_C.*degrat40yr_onshorewind(4,1).*utilize_ratio2060_2+cost_onshorewind_alone(:,9)./CO2_C.*degrat40yr_onshorewind(4,1)+ cost_onshorewind_alone(:,10)./CO2_C); % Mton CO2
CO2_offshorewind_all_utilize_trans_storage(:,1)=emissionfactor./fossilfuel_emissionfactor.*(cost_offshorewind_alone(:,5)./CO2_C.*degrat40yr_offshorewind(4,1).*utilize_ratio2060_2); % Mton CO2

CO2_all_utilize_trans_storage = CO2_PV_all_utilize_trans_storage; % Mton CO2
[m,n]=find(CO2_onshorewind_all_utilize_trans_storage~=0);
CO2_all_utilize_trans_storage(sub2ind(size(CO2_all_utilize_trans_storage), m, n))= CO2_onshorewind_all_utilize_trans_storage(sub2ind(size(CO2_onshorewind_all_utilize_trans_storage), m, n));
[m,n]=find(CO2_offshorewind_all_utilize_trans_storage~=0);
CO2_all_utilize_trans_storage(sub2ind(size(CO2_all_utilize_trans_storage), m, n))= CO2_offshorewind_all_utilize_trans_storage(sub2ind(size(CO2_offshorewind_all_utilize_trans_storage), m, n));

rmb2us=1/6.8967; % RMB to USD2019

%%
[m,n]=find(ixx_2==1);
EP(m,1) = 0.043037244; % coal,  $/kWh
[m,n]=find(ixx_2==2);
EP(m,1) = 0.058369861; % gas,  $/kWh
[m,n]=find(ixx_2==3);
EP(m,1) = 0.141408336; % oil,  $/kWh
CO2_all_utilize_trans_storage(CO2_all_utilize_trans_storage>0)=0;
B_utilize_trans_storage = (sum(lcoeeall_utilize_trans_storage(:,1:9),2)+Cost_storage+Cost_mechanical-phhall_utilize_trans_storage(:,1).*EP*10^3)./(-CO2_all_utilize_trans_storage); %USD/t CO2

%
IX=[1:1:numpowerunit]';
CO2_year_abated = CO2_PV_year_abated; % Mton CO2
[m,n]=find(CO2_onshorewind_year_abated~=0);
CO2_year_abated(sub2ind(size(CO2_year_abated), m, n))= CO2_onshorewind_year_abated(sub2ind(size(CO2_onshorewind_year_abated), m, n));
[m,n]=find(CO2_offshorewind_year_abated~=0);
CO2_year_abated(sub2ind(size(CO2_year_abated), m, n))= CO2_offshorewind_year_abated(sub2ind(size(CO2_offshorewind_year_abated), m, n));

CO2_year_ls = CO2_PV_year_ls; % Mton CO2
[m,n]=find(CO2_onshorewind_year_ls~=0);
CO2_year_ls(sub2ind(size(CO2_year_ls), m, n))= CO2_onshorewind_year_ls(sub2ind(size(CO2_onshorewind_year_ls), m, n));
[m,n]=find(CO2_offshorewind_year_ls~=0);
CO2_year_ls(sub2ind(size(CO2_year_ls), m, n))= CO2_offshorewind_year_ls(sub2ind(size(CO2_offshorewind_year_ls), m, n));
CO2_year_utilize_trans_storage = CO2_year_abated.*utilize_ratio2060_2+CO2_year_ls;
CO2_year_utilize_trans_storage = emissionfactor/fossilfuel_emissionfactor.*CO2_year_utilize_trans_storage;

CO2all_c_utilize_trans_storage = zeros(size(IX,1),1);
CO2all_c_utilize_trans_storage= CO2_year_utilize_trans_storage(sub2ind(size(CO2_year_utilize_trans_storage), IX, ones(size(IX,1),1)));
CO2all_c_utilize_trans_storage(CO2all_c_utilize_trans_storage>0)=0;
[B_utilize_trans_storage_IX,IX]=sort(B_utilize_trans_storage);
CP11122 = CP111(IX,1);
CP_cumsum_all=cumsum(CP11122);
[~,Index] = min(abs(CP_cumsum_all-8));

ID_Plant(:,1) =optpowerunit_IX(:,40);
ID_Plant(:,2) =optpowerunit_IX(:,35); % type ,1-3
ID_Plant_IX = ID_Plant(IX,:);
ID_Plant_IX_8TW = ID_Plant_IX(1:Index,:);
save('H:\China C neutrality\ANS\CO2_Battery.mat','CO2all_c_utilize_trans_storage'); 
save('H:\China C neutrality\ANS\B_Battery.mat','B_utilize_trans_storage'); 

for i = 1:size(B_utilize_trans_storage,1)
    CO2all_c_utilize_trans_storage_IX2(i,1) = CO2all_c_utilize_trans_storage(IX(i,1),1);
    B_utilize_trans_storage_IX2(i,1) = B_utilize_trans_storage(IX(i,1),1);
end

CO2all_c2_utilize_trans_storage_IX2=cumsum(CO2all_c_utilize_trans_storage_IX2);
CO2all2_c22_utilize_trans_storage_IX2(2,1) = 0;
CO2all2_c22_utilize_trans_storage_IX2(2,1) =  CO2all_c2_utilize_trans_storage_IX2(1);

B22_utilize_trans_storage_IX2(1,1) =  B_utilize_trans_storage_IX2(1);
B22_utilize_trans_storage_IX2(2,1) =  B_utilize_trans_storage_IX2(1);

for i = 2:size(CO2all_c2_utilize_trans_storage_IX2,1)
    CO2all2_c22_utilize_trans_storage_IX2(i*2-1,1) =  CO2all_c2_utilize_trans_storage_IX2(i-1,1);
    CO2all2_c22_utilize_trans_storage_IX2(i*2,1) =  CO2all_c2_utilize_trans_storage_IX2(i,1);
    B22_utilize_trans_storage_IX2(i*2-1,1) =  B_utilize_trans_storage_IX2(i,1);
    B22_utilize_trans_storage_IX2(i*2,1) =  B_utilize_trans_storage_IX2(i,1);
end  
CO2all2_c22_utilize_trans_storage_IX2_central = CO2all2_c22_utilize_trans_storage_IX2;
B22_utilize_trans_storage_IX2_central = B22_utilize_trans_storage_IX2;
save('H:\China C neutrality\ANS\CO2_Battery_ix.mat','CO2all2_c22_utilize_trans_storage_IX2'); 
save('H:\China C neutrality\ANS\B_Battery_ix.mat','B22_utilize_trans_storage_IX2'); 

