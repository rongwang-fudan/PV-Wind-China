tic
clear;
load('G:\China C neutrality\PV_power potential\ANS_PV1\tranmission_lines.dat','-mat');  % lines
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
lines(numlines+1:end,:)=[];

load('G:\China C neutrality\PV_power potential\ANS_PV1\optpowerunit_PV.mat'); % 
load('G:\China C neutrality\PV_power potential\ANS_PV1\powerunit_IX_PV.mat'); % 
optpowerunit_PV(:,35) = 1;
optpowerunit_PV(:,40) = powerunit_IX_PV; % power plant ID

load('G:\China C neutrality\Onshore wind_power potential\ANS_ONS1\optpowerunit_onshorewind.mat'); % 
load('G:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit_IX_onshorewind.mat'); % 

optpowerunit_onshorewind(:,35) = 2;
optpowerunit_onshorewind(:,40) = powerunit_IX_onshorewind; % power plant ID

load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\optpowerunit_offshorewind.mat'); % 
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\powerunit_IX_offshorewind.mat'); % 
optpowerunit_offshorewind(:,20)=optpowerunit_offshorewind(:,8);
optpowerunit_offshorewind(:,30)=optpowerunit_offshorewind(:,3)/1000; %MW
optpowerunit_offshorewind(:,35) = 3;
optpowerunit_offshorewind(:,40) = powerunit_IX_offshorewind; % power plant ID

optpowerunit = [optpowerunit_PV;optpowerunit_onshorewind;optpowerunit_offshorewind];
[B,IX]=sort(optpowerunit(:,20),1);
numpowerunit = size(optpowerunit,1);
for i=1:numpowerunit
    i2=IX(i);
    powerunit_IX(i,1)=i2;
    optpowerunit_IX(i,1:40)=optpowerunit(i2,1:40); % lat lon
end
% PV
LR_PV1 = 0.2; 
LR_PV2 = 0.18; 
CP0_PV = 253000; % cumulative capacity potential in China in 2020, MW
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
if optpowerunit_IX(1,35)==1
    cost_PV(1,11:19) =  optpowerunit_IX(1,11:19);
end
for i = 2: numpowerunit
    if optpowerunit_IX(i,35)==1
        cost_PV(i,11:19) =  cost_PV(i-1,11:19) + optpowerunit_IX(i,11:19);
    else
        cost_PV(i,11:19) =  cost_PV(i-1,11:19);
    end    
end


% onshorewind
LR_onshorewind1= 0.074; 
LR_onshorewind2= 0.18;
CP0_onshorewind = 272010; % MW
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
if optpowerunit_IX(1,35)==2
    cost_onshorewind(1,11:19) =  optpowerunit_IX(1,11:19);
end
for i = 2: numpowerunit
    if optpowerunit_IX(i,35)==2
        cost_onshorewind(i,11:19) =  cost_onshorewind(i-1,11:19) + optpowerunit_IX(i,11:19);
    else
        cost_onshorewind(i,11:19) =  cost_onshorewind(i-1,11:19);
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
cost_storage_offshorewind = zeros(numpowerunit,1);
if optpowerunit_IX(1,35)==3
    cost_offshorewind(1,6:10) =  optpowerunit_IX(1,6:10);
end
for i = 2: numpowerunit
    if optpowerunit_IX(i,35)==3
        cost_offshorewind(i,6:10) =  cost_offshorewind(i-1,6:10) + optpowerunit_IX(i,6:10);
    else
        cost_offshorewind(i,6:10) =  cost_offshorewind(i-1,6:10);
    end    
end
%%
load('G:\China C neutrality\Data\powerdemand_weekday_CN2060all_pro.mat')  % TWh/h  
load('G:\China C neutrality\Data\powerdemand_weekday_CN2060_pro')% TWh/h  
rr = 0.58/(0.42*0.342862283+0.58);
powerdemand_monhour2060all = powerdemand_weekday_CN2060all.*rr;
p_all2060=sum(sum(powerdemand_monhour2060all))/168*8760; % 	TWh/year

P_PV111= zeros(numpowerunit,1);
P_onshorewind111= zeros(numpowerunit,1);
P_offshorewind111= zeros(numpowerunit,1);
for i = 1: numpowerunit
    if optpowerunit_IX(i,35)==1
        P_PV111(i,1) =  optpowerunit_IX(i,1);
    else if optpowerunit_IX(i,35)==2
            P_onshorewind111(i,1) =  optpowerunit_IX(i,1);
        else if optpowerunit_IX(i,35)==3
                P_offshorewind111(i,1) =  optpowerunit_IX(i,1);
            end
        end
    end
end
P_ratio_3type(:,1) = cumsum(P_PV111)./p_all2060;
P_ratio_3type(:,2) = cumsum(P_onshorewind111)./p_all2060;
P_ratio_3type(:,3) = cumsum(P_offshorewind111)./p_all2060;

(sum(P_PV111)+sum(P_onshorewind111)+sum(P_offshorewind111))/p_all2060
(sum(P_PV111)+sum(P_onshorewind111)+sum(P_offshorewind111))/(sum(sum(powerdemand_weekday_CN2060))/168*8760/(1+0.02)^41)

optpowerunit1= zeros(numpowerunit,30);
optpowerunit1(1,30)=optpowerunit_IX(1,30);
optpowerunit1(1,1)=optpowerunit_IX(1,1);

for i=2:numpowerunit
    optpowerunit1(i,30)=optpowerunit1(i-1,30)+optpowerunit_IX(i,30); % cumulative capacity
    optpowerunit1(i,1)=optpowerunit1(i-1,1)+optpowerunit_IX(i,1); % electricity generation, TWh/year
end
optpowerunit_IX(:,36)=optpowerunit1(:,30)/optpowerunit1(end,30);
pe_ratio = optpowerunit1(:,1)/p_all2060;


costmin=1e18;
ratio_onshorewind = 0.704960835509138;
ratio_PV = 0.497159090909091;
module_price_offshore = 1.008; % USD/W

for p1=0:1:100
    display(p1);
    pene=zeros(4,2);
    pene(1,1)=p1/100; % peneratration ratio in 2020-2030
    for p2=0:1:(100-p1)
        pene(2,1)=(p1+p2)/100; % peneratration ratio in 2030-2040        
        for p3=0:1:(100-p1-p2)
            pene(3,1)=(p1+p2+p3)/100; % peneratration ratio in 2040-2050            
            pene(4,1)=1; % peneratration ratio in 2050-2060
                costall_PV=0;
                costall_onshorewind=0;
                costall_offshorewind=0;
                cumulativecapacity=0;
                unit1=1;
                for yy=1:4
                    idx=find(optpowerunit_IX(:,36)>=pene(yy,1)); % percentage of energy
                    unit2(yy,1)=idx(1);
                    if unit2(yy,1)>unit1
                        if unit1==1
                            % PV
                            costall_PV=costall_PV+cost_PV(unit2(yy,1),14) *ratio_PV +cost_PV(unit2(yy,1),14) *(1-ratio_PV)+cost_PV(unit2(yy,1),16) +sum(cost_PV(unit2(yy,1),11:13),2)+ cost_PV(unit2(yy,1),15); % cost million USD 
                            
                            % onshorewind
                            costall_onshorewind=costall_onshorewind+cost_onshorewind(unit2(yy,1),14) *ratio_onshorewind+cost_onshorewind(unit2(yy,1),16) +cost_onshorewind(unit2(yy,1),14) *(1-ratio_onshorewind)+sum(cost_onshorewind(unit2(yy,1),11:13),2); % cost million USD 

                            % offshorewind
                            costall_offshorewind=costall_offshorewind+cost_offshorewind(unit2(yy,1),6)+cost_offshorewind(unit2(yy,1),7); % cost million USD 
                        else
                            % PV
                            costall_PV=costall_PV+(cost_PV(unit2(yy,1),14)-cost_PV(unit1-1,14))*learnprice_PV(unit1)+(cost_PV(unit2(yy,1),16)-cost_PV(unit1-1,16))+(sum(cost_PV(unit2(yy,1),11:13),2)-sum(cost_PV(unit1-1,11:13),2))*learnprice_PV_2(unit1)+(cost_PV(unit2(yy,1),15)-cost_PV(unit1-1,15))*learnprice_PV(unit1); % cost million USD 

                            % onshorewind
                            costall_onshorewind=costall_onshorewind+(cost_onshorewind(unit2(yy,1),14)-cost_onshorewind(unit1-1,14)) *learnprice_onshorewind(unit1)+(cost_onshorewind(unit2(yy,1),16)-cost_onshorewind(unit1-1,16))+(sum(cost_onshorewind(unit2(yy,1),11:13),2)-sum(cost_onshorewind(unit1-1,11:13),2))*learnprice_onshorewind_2(unit1); % cost million USD 

                            % offshorewind
                            % 1  module cost
                            aa= module_price_offshore.*(CP0_offshorewind_cum(unit2(yy,1),1)-CP0_offshorewind_cum(unit1-1,1))*learnprice_offshorewind(unit1); % cost million USD 
                            % 2  others
                            aa2= ((cost_offshorewind(unit2(yy,1),6)-cost_offshorewind(unit1-1,6))- module_price_offshore.*(CP0_offshorewind_cum(unit2(yy,1),1)-CP0_offshorewind_cum(unit1-1,1)))*learnprice_offshorewind_2(unit1); % cost million USD 
                            % 9  O&M cost
                            aa3 = (cost_offshorewind(unit2(yy,1),7)-cost_offshorewind(unit1-1,7)); % cost million USD 
                            costall_offshorewind=costall_offshorewind+aa+aa2+aa3; % cost million USD 
                        end
                        pene(yy,2)=optpowerunit_IX(unit2(yy,1),36); % peneratration ratio
                        pene_real(yy,2)=pe_ratio(unit2(yy,1),1); % The ratio of electricity generation to total electricity demand in 2060
                        unit1=unit2(yy,1)+1;
                    end
                end
                costall = costall_PV + costall_onshorewind + costall_offshorewind;
                if costall<costmin
                    costmin=costall;
                    unitmin=unit2;
                end
        end
    end
end
save('G:\China C neutrality\Data\ANS\unitmin_module.mat','unitmin');

