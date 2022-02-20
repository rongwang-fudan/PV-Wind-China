tic
clear;
load('G:\China C neutrality\PV_power potential\ANS_PV1\tranmission_lines.dat','-mat');  % lines
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
lines(numlines+1:end,:)=[];
load('G:\China C neutrality\Data\ANS\unitmin_module.mat'); % The number of PV,onshore wind and offshore wind power plants at each stage of construction at minimum cost

load('G:\China C neutrality\Data\ANS\utilize_ratio2060_plant_alone.mat') % power use efficiency without inter-regional electric transport and storage
utilize_ratio2060_plant_alone1=utilize_ratio2060_plant_alone;
load('G:\China C neutrality\Data\ANS\utilize_ratio2060_trans_plant_alone.mat')  % power use efficiency with inter-regional electric transport but without storage
utilize_ratio2060_trans_plant_alone1=utilize_ratio2060_trans_plant_alone;
load('G:\China C neutrality\Data\ANS\utilize_ratio2060_trans_storage_plant_alone.mat')  % power use efficiency with inter-regional electric transport and storage
utilize_ratio2060_trans_storage_plant_alone1=utilize_ratio2060_trans_storage_plant_alone;
load('G:\China C neutrality\Data\ANS\storage_max_plant.mat')   % Maximum charging power, TWh/h  
load('G:\China C neutrality\Data\ANS\storage_year_plant.mat')  % Average annual electricity storage, TWh/year
load('G:\China C neutrality\Data\ANS\cost_trans_IX2.mat')% transmission cost

lifetime_s = 5; % year
Ip = storage_max_plant*10^9; % kW  power capacity
Ie = storage_year_plant*10^9*lifetime_s/8000; % kWh  energy capacity

[m,n] =find(Ie==0);
Ie_cum = cumsum(Ie)/10^6;% GWh
Index_Ie(1,1) = m(end)+1; % The position to start the storage

for i = 1:floor(Ie_cum(end)/50)
    [m,n]=find(Ie_cum/50>i);
    Index_Ie(i+1,1) = m(1);
end
Index_Ie(floor(Ie_cum(end)/50)+1,1) = 3845;

p_dis = storage_year_plant*10^9*0.99*0.85; % annual discharge, kWh/year 
p_char = storage_year_plant*10^9; % annual charge, kWh/year 
Co = 1.5/1000; % $/kWh
Cp=595.73; % $/kW
discount=0.07; % per year
discount1yr=0;
for t=1:40
    discount1yr=discount1yr+1/(1+discount)^(t);
end
%
Ce = 357; % 2018$/kWh
LR_sto1 = 0; 
LR_sto2 = 0.3; 
CP0_sto = sum(sum(Ie))/10; % kW
CP_sto = zeros(size(Ie,1),1);
CP_sto(1,1) = 2.91*10^6; %kw
for i = 1: size(Ie,1)
        CP_sto(i+1,1) =  CP_sto(i,1) + Ie(i,1);
end
r_sto = CP_sto./CP0_sto;
for i = 1:size(r_sto,1)
    if CP_sto(i,1)<=CP0_sto
        LR_sto(i,1) = LR_sto1;
        learnprice_sto(i,1) = (CP_sto(i,1)./CP0_sto).^(log(1- LR_sto(i,1))/log(2));
    else
        LR_PV_2(i,1) = LR_sto2;
        learnprice_sto(i,1) = (CP_sto(i,1)./CP0_sto).^(log(1- LR_PV_2(i,1))/log(2));      
    end
end

Cost_storage1 = (Ce.*Ie+Cp.*Ip+Co*(p_dis+p_char)*discount1yr)/10^6.*learnprice_sto(3844,1)*40/lifetime_s*discount1yr/40; %million $
Cost_ele1=sum(sum(Cost_storage1))./(sum(sum(Ie)))*10^6*50;% million $/50GWh    
clear discount1yr
%%
CO2_C=0.2727;
lifetime_power=40;
lifetime_inverter=10; % renewed per 10 years
discount=0.07; % per year
OMratio_PV = 1/lifetime_power*0.4;
OMratio_wind=1/lifetime_power*1.2;
discountinverter=0;
for t=1:floor(lifetime_power/lifetime_inverter)
    discountinverter=discountinverter+1/(1+discount)^((t-1)*lifetime_inverter);
end

discount40yr=zeros(4,1);
for t=1:lifetime_power
    discount40yr(4,1)=discount40yr(4,1)+1/(1+discount)^(t-1);
end
for t=1:30
    discount40yr(3,1)=discount40yr(3,1)+1/(1+discount)^(t-1);
end
for t=1:20
    discount40yr(2,1)=discount40yr(2,1)+1/(1+discount)^(t-1);
end
for t=1:10
    discount40yr(1,1)=discount40yr(1,1)+1/(1+discount)^(t-1);
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
for t=1:30
    degrat40yr_PV(3,1)=degrat40yr_PV(3,1)+(1-degration_PV)^(t-1)/(1+discount)^(t-1);
    degrat40yr_onshorewind(3,1)=degrat40yr_onshorewind(3,1)+(1-degration_onshorewind)^(t-1)/(1+discount)^(t-1);
    degrat40yr_offshorewind(3,1)=degrat40yr_offshorewind(3,1)+(1-degration_offshorewind)^(t-1)/(1+discount)^(t-1);
end
for t=1:20
    degrat40yr_PV(2,1)=degrat40yr_PV(2,1)+(1-degration_PV)^(t-1)/(1+discount)^(t-1);
    degrat40yr_onshorewind(2,1)=degrat40yr_onshorewind(2,1)+(1-degration_onshorewind)^(t-1)/(1+discount)^(t-1);
    degrat40yr_offshorewind(2,1)=degrat40yr_offshorewind(2,1)+(1-degration_offshorewind)^(t-1)/(1+discount)^(t-1);
end
for t=1:10
    degrat40yr_PV(1,1)=degrat40yr_PV(1,1)+(1-degration_PV)^(t-1)/(1+discount)^(t-1);
    degrat40yr_onshorewind(1,1)=degrat40yr_onshorewind(1,1)+(1-degration_onshorewind)^(t-1)/(1+discount)^(t-1);
    degrat40yr_offshorewind(1,1)=degrat40yr_offshorewind(1,1)+(1-degration_offshorewind)^(t-1)/(1+discount)^(t-1);
end

load('G:\China C neutrality\PV_power potential\ANS_PV1\optpowerunit_PV.mat'); % 
load('G:\China C neutrality\PV_power potential\ANS_PV1\powerunit_IX_PV.mat'); % 
load('G:\China C neutrality\PV_power potential\ANS_PV1\tranmission_lines_IX.mat');  % lines_IX
lines_IX(size(optpowerunit_PV,1)+1:end,:)=[];
lines_IX_PV = lines_IX;
load('G:\China C neutrality\PV_power potential\ANS_PV1\powerunit_num_IX_PV.mat');  % 
load('G:\China C neutrality\PV_power potential\ANS_PV1\unitid_lcoe.dat','-mat'); 
unitid_lcoe_PV = unitid_lcoe;
optpowerunit_PV(:,35) = 1;
optpowerunit_PV(:,40) = powerunit_IX_PV; 

load('G:\China C neutrality\Onshore wind_power potential\ANS_ONS1\optpowerunit_onshorewind.mat'); % 
load('G:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit_IX_onshorewind.mat'); % 
load('G:\China C neutrality\Onshore wind_power potential\ANS_ONS1\tranmission_lines_IX.mat');  % lines_IX
lines_IX(size(optpowerunit_onshorewind,1)+1:end,:)=[];
lines_IX_onshorewind = lines_IX;
load('G:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit_num_IX_onshorewind.mat');  % 
optpowerunit_onshorewind(:,35) = 2;
optpowerunit_onshorewind(:,40) = powerunit_IX_onshorewind;

load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\optpowerunit_offshorewind.mat'); % 
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\powerunit_IX_offshorewind.mat'); % 
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\tranmission_lines_IX.mat');  % lines_IX
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
CP_cumsum_all = cumsum(optpowerunit_IX(:,30))/10^6;
CF=sum(optpowerunit_IX(:,1))./sum(optpowerunit_IX(:,30))/8760*10^6;
[B,IX]=sort(optpowerunit_IX(:,20),1);

%%
for i = 1:floor(Ie_cum(end)/50)
    r1=sum(optpowerunit_IX(Index_Ie(i,1):Index_Ie(i+1,1)-1,1).*utilize_ratio2060_trans_storage_plant_alone(Index_Ie(i,1):Index_Ie(i+1,1)-1,1))/sum(optpowerunit_IX(Index_Ie(i,1):Index_Ie(i+1,1)-1,1));
    utilize_ratio2060_trans_storage_plant_alone(Index_Ie(i,1):Index_Ie(i+1,1)-1,1)=r1;
    Cost_storage(Index_Ie(i,1):Index_Ie(i+1,1)-1,1) = (optpowerunit_IX(Index_Ie(i,1):Index_Ie(i+1,1)-1,1).* optpowerunit_IX(Index_Ie(i,1):Index_Ie(i+1,1)-1,20)*1000.*degrat40yr_offshorewind(4,1)).*Cost_ele1./sum(optpowerunit_IX(Index_Ie(i,1):Index_Ie(i+1,1)-1,1).* optpowerunit_IX(Index_Ie(i,1):Index_Ie(i+1,1)-1,20)*1000.*degrat40yr_offshorewind(4,1)); %million
end
clear Cost_storage1

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
CO2_PV_alone = zeros(numpowerunit,3);
if optpowerunit_IX(1,35)==1
    cost_PV(1,11:19) =  optpowerunit_IX(1,11:19);
    cost_PV_alone(1,1:30) = optpowerunit_IX(1,1:30);
    CO2_PV_alone(1,1:3)  = optpowerunit_IX(1,8:10);
end
for i = 2: numpowerunit
    if optpowerunit_IX(i,35)==1
        cost_PV(i,11:19) =  cost_PV(i-1,11:19) + optpowerunit_IX(i,11:19);
        cost_PV_alone(i,1:30) =  optpowerunit_IX(i,1:30);
        CO2_PV_alone(i,1:3)  = optpowerunit_IX(i,8:10);
    else
        cost_PV(i,11:19) =  cost_PV(i-1,11:19);
        cost_PV_alone(i,1:30) =  0;
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
CO2_onshorewind_alone = zeros(numpowerunit,3);
if optpowerunit_IX(1,35)==2
    cost_onshorewind(1,11:19) =  optpowerunit_IX(1,11:19);
    cost_onshorewind_alone(1,1:30) =  optpowerunit_IX(1,1:30);
    CO2_onshorewind_alone(1,1:3) =  optpowerunit_IX(1,8:10);
end
for i = 2: numpowerunit
    if optpowerunit_IX(i,35)==2
        cost_onshorewind(i,11:19) =  cost_onshorewind(i-1,11:19) + optpowerunit_IX(i,11:19);
        cost_onshorewind_alone(i,1:30) =  optpowerunit_IX(i,1:30);
        CO2_onshorewind_alone(i,1:3) =  optpowerunit_IX(i,8:10);
    else
        cost_onshorewind(i,11:19) =  cost_onshorewind(i-1,11:19);
        cost_onshorewind_alone(i,1:30) =  0;
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
if optpowerunit_IX(1,35)==3
    cost_offshorewind(1,6:10) =  optpowerunit_IX(1,6:10);
    cost_offshorewind_alone(1,1:30) =  optpowerunit_IX(1,1:30);
end
for i = 2: numpowerunit
    if optpowerunit_IX(i,35)==3
        cost_offshorewind(i,6:10) =  cost_offshorewind(i-1,6:10) + optpowerunit_IX(i,6:10);
        cost_offshorewind_alone(i,1:30) =  optpowerunit_IX(i,1:30);        
    else
        cost_offshorewind(i,6:10) =  cost_offshorewind(i-1,6:10);
        cost_offshorewind_alone(i,1:30) =  0;
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

sto_PV(:,1) = cost_PV_alone(:,1).*(1-utilize_ratio2060_trans_plant_alone); 
phh_PV_alone(:,1) = cost_PV_alone(:,1); 
sto_PV_sto2(:,1) = cost_PV_alone(:,1).*(1-utilize_ratio2060_trans_storage_plant_alone); 

% power generation
phh_offshorewind(:,1) = cost_offshorewind_alone(:,1).*degrat40yr_offshorewind(4,1); 

% power generation
phh_onshorewind_utilize_trans_storage(:,1) = cost_onshorewind_alone(:,1).*degrat40yr_onshorewind(4,1).*utilize_ratio2060_trans_storage_plant_alone; 

module_price_offshore = 1.008; % USD/W
% PV
% 1  fixed system cost
lcoee_PV(:,1) = cost_PV_alone(:,14)./disrate_PV(4,1) .*learnprice_PV(3844,1) + cost_PV_alone(:,15)/disrate_inverter_PV(4,1).*discountinverter .*learnprice_PV(3844,1) ; % cost million USD 
% 2  substation
lcoee_PV(:,2) = cost_PV_alone(:,12)./disrate_PV(4,1).*learnprice_PV_2(3844,1) ; % cost million USD 
% 3  Power plant cable
lcoee_PV(:,3) = cost_PV_alone(:,13)./disrate_PV(4,1).*learnprice_PV_2(3844,1); % cost million USD 
% 5  connection to national grid
lcoee_PV(:,5) = cost_PV_alone(:,11)./disrate_PV(4,1).*learnprice_PV_2(3844,1); % cost million USD 
% 9  O&M cost
lcoee_PV(:,9) = lcoee_PV(:,4).*(disrate_PV(4,1)-1)+(sum(lcoee_PV(:,1:3),2)+lcoee_PV(:,5)).*(disrate_PV(4,1)-1); % cost million USD 

% onshorewind
% 1  fixed system cost
lcoee_onshorewind(:,1) = cost_onshorewind_alone(:,14)./disrate_wind(4,1).*learnprice_onshorewind(3844,1); % cost million USD 
% 2  substation
lcoee_onshorewind(:,2) = cost_onshorewind_alone(:,12)./disrate_wind(4,1).*learnprice_onshorewind_2(3844,1); % cost million USD 
% 3  Power plant cable
lcoee_onshorewind(:,3) = cost_onshorewind_alone(:,13)./disrate_wind(4,1).*learnprice_onshorewind_2(3844,1); % cost million USD 
% 5  connection to national grid
lcoee_onshorewind(:,5) = cost_onshorewind_alone(:,11)./disrate_wind(4,1).*learnprice_onshorewind_2(3844,1); % cost million USD 
% 9  O&M cost
lcoee_onshorewind(:,9) = lcoee_onshorewind(:,4).*(disrate_wind(4,1)-1)+(sum(lcoee_onshorewind(:,1:3),2)+lcoee_onshorewind(:,5)).*(disrate_wind(4,1)-1); % cost million USD 

% offshorewind
% 1  module cost
lcoee_offshorewind(:,1) = module_price_offshore.*CP_offshorewind_alone(:,1).*learnprice_offshorewind(3844,1); % cost million USD 
% 2  others
lcoee_offshorewind(:,2) = (cost_offshorewind_alone(:,6)./disrate_wind(4,1)-module_price_offshore.*CP_offshorewind_alone(:,1)).*learnprice_offshorewind_2(3844,1) ; % cost million USD 
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

CO2_PV_all_utilize_trans_storage(:,1)= cost_PV_alone(:,8)./CO2_C.*degrat40yr_PV(4,1).*utilize_ratio2060_trans_storage_plant_alone+ cost_PV_alone(:,9)./CO2_C.*degrat40yr_PV(4,1)+ cost_PV_alone(:,10)./CO2_C; % Mton CO2
CO2_onshorewind_all_utilize_trans_storage(:,1)=cost_onshorewind_alone(:,8)./CO2_C.*degrat40yr_onshorewind(4,1).*utilize_ratio2060_trans_storage_plant_alone+cost_onshorewind_alone(:,9)./CO2_C.*degrat40yr_onshorewind(4,1)+ cost_onshorewind_alone(:,10)./CO2_C; % Mton CO2
CO2_offshorewind_all_utilize_trans_storage(:,1)=cost_offshorewind_alone(:,5)./CO2_C.*degrat40yr_offshorewind(4,1).*utilize_ratio2060_trans_storage_plant_alone; % Mton CO2

CO2_all_utilize_trans_storage = CO2_PV_all_utilize_trans_storage; % Mton CO2
[m,n]=find(CO2_onshorewind_all_utilize_trans_storage~=0);
CO2_all_utilize_trans_storage(sub2ind(size(CO2_all_utilize_trans_storage), m, n))= CO2_onshorewind_all_utilize_trans_storage(sub2ind(size(CO2_onshorewind_all_utilize_trans_storage), m, n));
[m,n]=find(CO2_offshorewind_all_utilize_trans_storage~=0);
CO2_all_utilize_trans_storage(sub2ind(size(CO2_all_utilize_trans_storage), m, n))= CO2_offshorewind_all_utilize_trans_storage(sub2ind(size(CO2_offshorewind_all_utilize_trans_storage), m, n));

%%
rmb2us=1/6.8967; % RMB to USD2019
EP = 0.15*rmb2us; %$/kWh
B_utilize_trans_storage = (sum(lcoeeall_utilize_trans_storage(:,1:9),2)+Cost_storage-phhall_utilize_trans_storage(:,1)*EP*10^3)./(-CO2_all_utilize_trans_storage); %USD/t CO2

%%
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

CO2_year_utilize_trans_storage = CO2_year_abated.*utilize_ratio2060_trans_storage_plant_alone+CO2_year_ls;

%%
CO2all_c_utilize_trans_storage = zeros(size(IX,1),1);
CO2all_c_utilize_trans_storage= CO2_year_utilize_trans_storage(sub2ind(size(CO2_year_utilize_trans_storage), IX, ones(size(IX,1),1)));
[B_utilize_trans_storage_IX,IX]=sort(B_utilize_trans_storage);
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

%
figure
plot(-CO2all2_c22_utilize_trans_storage_IX2./10^3,B22_utilize_trans_storage_IX2) 
hold on
% axis([0 max(max(-CO2all2_c22_utilize_trans_storage_IX2./10^3)) -20 220]);
xlabel('Emissions reduced by transiting fossil fuel plants to PV and wind power plants in 2060 (Gt CO2 y-1)')
ylabel('MAC (2019 $ (t CO2)-1)')

