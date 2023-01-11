tic
clear

load('H:\China C neutrality\ANS\CO2_mecha_xz.mat'); % USD2019/kWh
load('H:\China C neutrality\ANS\B_mecha_xz.mat'); % USD2019/kWh
load('H:\China C neutrality\ANS\LCOEE_mecha_xz.mat'); % USD2019/kWh
CO2_mecha = CO2all_c_utilize_trans_storage;
B_mecha = B_utilize_trans_storage;
LCOE_mecha = LCOEE_all_utilize_trans_storage;

load('H:\China C neutrality\ANS\CO2_Battery_xz.mat'); % USD2019/kWh
load('H:\China C neutrality\ANS\B_Battery_xz.mat'); % USD2019/kWh
load('H:\China C neutrality\ANS\LCOEE_Battery_xz.mat'); % USD2019/kWh
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

load('H:\China C neutrality\ANS\storage_max_plant_xz.mat')   % Maximum charging power, TWh/h  
load('H:\China C neutrality\ANS\storage_year_plant_xz.mat')  % Average annual electricity storage, TWh/year
Ip = storage_max_plant/(0.99*0.85*0.983*0.967)*10^9; % kW  power capacity
[m,n]=find(cho==1);
sum(Ip(m))
[m,n]=find(cho==-1);
sum(Ip(m))

[B_IX,IX]=sort(B);
CO2_IX = CO2(IX,1);
[B_Battery_IX,IX1]=sort(B_Battery);
CO2_Battery_IX = CO2_Battery(IX1,1);
[B_mecha_IX,IX2]=sort(B_mecha);
CO2_mecha_IX = CO2_mecha(IX2,1);

save('H:\China C neutrality\ANS\cho_xz.mat','cho'); % 
