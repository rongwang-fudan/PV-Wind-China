tic
clear;
rmb2us=1/6.8967; % RMB to USD2019
fossilfuel_emissionfactor=0.827; % tCO2/MWh 
CO2_C=0.2727;
lifetime_power=40;
discount=0.07; % per year
discount1yr=0;
for t=1:lifetime_power
    discount1yr=discount1yr+1/(1+discount)^(t-1);
end
% PV:
OMratio_majorline_PV=1/lifetime_power*0.2;
OMratio_substation_PV=1/lifetime_power*0.2;
% wind
OMratio_majorline_wind=1/lifetime_power*1.2;
OMratio_substation_wind=1/lifetime_power*1.2;

load('G:\China C neutrality\PV_power potential\ANS_PV1\tranmission_lines.dat','-mat');  % lines
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
lines(numlines+1:end,:)=[];

load('G:\China C neutrality\PV_power potential\ANS_PV1\optpowerunit_PV.mat'); % 
load('G:\China C neutrality\PV_power potential\ANS_PV1\powerunit_IX_PV.mat'); % 
load('G:\China C neutrality\PV_power potential\ANS_PV1\tranmission_lines_IX.mat');  % lines_IX
lines_IX(size(optpowerunit_PV,1)+1:end,:)=[];
lines_IX_PV = lines_IX;
load('G:\China C neutrality\PV_power potential\ANS_PV1\powerunit_num_IX_PV.mat');  % 
optpowerunit_PV(:,35) = 1;
optpowerunit_PV(:,40) = powerunit_IX_PV; % power plant type ID

load('G:\China C neutrality\Onshore wind_power potential\ANS_ONS1\optpowerunit_onshorewind.mat'); % 
load('G:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit_IX_onshorewind.mat'); % 
load('G:\China C neutrality\Onshore wind_power potential\ANS_ONS1\tranmission_lines_IX.mat');  % lines_IX
lines_IX(size(optpowerunit_onshorewind,1)+1:end,:)=[];
lines_IX_onshorewind = lines_IX;
load('G:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit_num_IX_onshorewind.mat');  % 
optpowerunit_onshorewind(:,35) = 2;
optpowerunit_onshorewind(:,40) = powerunit_IX_onshorewind; % power plant type ID

load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\optpowerunit_offshorewind.mat'); % 
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\powerunit_IX_offshorewind.mat'); % 
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\tranmission_lines_IX.mat');  % lines_IX
lines_IX(size(powerunit_IX_offshorewind,1)+1:end,:)=[];
lines_IX_offshorewind = lines_IX;
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\powerunit_num_IX_offshorewind.mat');  % 
powerunit_num_IX_offshorewind(:,4) = 1;

optpowerunit_offshorewind(:,20)=optpowerunit_offshorewind(:,8);
optpowerunit_offshorewind(:,30)=optpowerunit_offshorewind(:,3)/1000; %MW
optpowerunit_offshorewind(:,35) = 3;
optpowerunit_offshorewind(:,40) = powerunit_IX_offshorewind; % power plant type ID

optpowerunit = [optpowerunit_PV;optpowerunit_onshorewind;optpowerunit_offshorewind];
lines_IX = [lines_IX_PV(:,1:15);lines_IX_onshorewind;lines_IX_offshorewind];
powerunit_num_IX = [powerunit_num_IX_PV;powerunit_num_IX_onshorewind;powerunit_num_IX_offshorewind];
[B,IX]=sort(optpowerunit(:,20),1);
numpowerunit = size(optpowerunit,1);
for i=1:numpowerunit
    i2=IX(i);
    powerunit_IX(i,1)=i2;
    optpowerunit_IX(i,1:40)=optpowerunit(i2,1:40); %  optpowerunit_IX(i2,1); electricity used by the county TWh / year
    lines_IX_IX(i,1:15)=lines_IX(i2,1:15); % lat lon
    powerunit_num_IX_IX(i,1:4)=powerunit_num_IX(i2,1:4); % lat lon
end
[B,IX]=sort(optpowerunit_IX(:,20),1);

line_IX_all = [lines(:,1:15);lines_IX_IX];
[m,n] = find(line_IX_all(:,9)==0);
line_IX_all(m,9) = 1362; % offshorewind, Shenzhen

%%
rate_energy=0.02; % per year 
pal = (7.6*0.599315068+2.5*0.922023182+1.8*0.7+2*0.6)/10*8760; % TWh/year
P_nuclear_utilize_trans_storage = 2.5*0.922023182/10*8760; %  TWh/year
P_hydro_utilize_trans_storage = 7.6*0.599315068/10*8760; %  TWh/year
P_BECCS_utilize_trans_storage = 1.8*0.7/10*8760; %  TWh/year
P_hydrogen_utilize_trans_storage = 2*0.6/10*8760; %  TWh/year
P_other_utilize_trans_storage = P_nuclear_utilize_trans_storage + P_hydro_utilize_trans_storage+P_BECCS_utilize_trans_storage+P_hydrogen_utilize_trans_storage; %   TWh/year

%
load('G:\China C neutrality\Data\P2019_others.mat')  % *10^8 kWh   bioenergy, nuclear, hydro, hydrogen energy
rr=P2019_others./sum(P2019_others);
P_BECCS2060 = P_BECCS_utilize_trans_storage.*rr(:,1);
P_nuclearS2060 = P_nuclear_utilize_trans_storage.*rr(:,2);
P_hydroS2060 = P_hydro_utilize_trans_storage.*rr(:,2);
P_hydrogenS2060 = P_hydrogen_utilize_trans_storage.*ones(34,1)*1/31;

load('G:\China C neutrality\Data\powerdemand_monhour2060all.mat') % TWh/h  
load('G:\China C neutrality\Data\powerdemand_monhour2060.mat')  % TWh/h  
rr = 0.58/(0.42*0.342862283+0.58);
powerdemand_monhour2060all = powerdemand_monhour2060all.*rr;
powerdemand_monhour = powerdemand_monhour2060all-(P_BECCS2060+P_nuclearS2060+P_hydroS2060+P_hydrogenS2060)'/8760;

regions = [2,3,2,5,6,6,7,6,3,1,6,4,1,1,4,2,2,4,0,3,5,5,7,5,2,2,3,3,6,5,7,7,2,7];
powerdemand_monhour_reg = zeros(288,7);
for i = 1:34
    if regions(i)~=0
        powerdemand_monhour_reg(:,regions(i)) = powerdemand_monhour(:,i) + powerdemand_monhour_reg(:,regions(i));
    end
end
powerdemand_monhour_CN = sum(powerdemand_monhour_reg,2);
rate_energy=0.02; % per year
powerdemand_monhour_CN2060 = sum(powerdemand_monhour_reg,2); 

load('G:\China C neutrality\Data\powergenerat_monhour_pv0928.mat')  % TWh/h  
load('G:\China C neutrality\Data\powergenerat_monhour_onshorewind0928.mat')   % TWh/h  
load('G:\China C neutrality\Data\powergenerat_monhour_offshorewind0928.mat')  % TWh/h  

%%
optpowerunit_IX_IX = optpowerunit_IX;
clear optpowerunit_IX
CP_cumsum_all = cumsum(optpowerunit_IX_IX(:,30))/10^6;
AA(:,1)=CP_cumsum_all;%TW
powergenerate_monhour_CN = zeros(288,1);
powergenerate_monhour_CN2 = zeros(288,1);
powergenerate_monhour_CN2_2060 = zeros(288,1);
storage = zeros(288,3844);
powergenerate_aaa=0;
aa=0;
aa_trans=0;
aa_trans_storage=0;
bb=0;
aa1=0;
bb1=0;
powergenerate_aaa = zeros(288,1);
powergenerate_monhour_CN2_2060_trans_storage = zeros(288,1);
bbbb = zeros(288,1);
aaaa = zeros(288,1);
aaaa_reg = zeros(288,7);
aaaa_trans = zeros(288,1);
aaaa_reg_trans = zeros(288,7);
aaaa_trans_storage = zeros(288,1);
aaaa_reg_trans_storage = zeros(288,7);
powergenerate_monhour_reg = zeros(288,7);
powergenerate_monhour_reg2_2060 = zeros(288,7);
powergenerate_monhour_reg2_2060_trans = zeros(288,7);
powergenerate_monhour_reg2_2060_trans_storage = zeros(288,7);
for i = 1:numpowerunit
    i2=i;
    type(i,1) = optpowerunit_IX_IX(i2,35);
    cy=line_IX_all(numlines+i2,9); % county of power unit
    reg=line_IX_all(numlines+i2,11); % region of power unit 1-7
    i3 =  optpowerunit_IX_IX(i2,31);
    if type(i,1)==1
        powergenerat_monhour(:,i) = powergenerat_monhour_pv(:,i3);
    else if type(i,1)==2
            powergenerat_monhour(:,i) = powergenerat_monhour_onshorewind(:,i3);
        else if type(i,1)==3
                powergenerat_monhour(:,i) = powergenerat_monhour_offshorewind(:,i3);
            end
        end
    end
    rrrr(i,1) = optpowerunit_IX_IX(i2,1)/(sum(powergenerat_monhour(:,i))/288*8760);
    powergenerat_monhour(:,i) = powergenerat_monhour(:,i).*optpowerunit_IX_IX(i2,1)/(sum(powergenerat_monhour(:,i))/288*8760);
    qwe(i,1)=sum(powergenerat_monhour(:,i))/288*8760 - optpowerunit_IX_IX(i2,1);
    powergenerat_monhour(:,i) = powergenerat_monhour(:,i).*optpowerunit_IX_IX(i,1)/(sum(powergenerat_monhour(:,i))'/288*8760);
    powergenerate_monhour_CN_22 = powergenerate_monhour_CN2;
    powergenerate_monhour_CN(:,1) = powergenerat_monhour(:,i) + powergenerate_monhour_CN(:,1);% TWh/h
    power_all(:,i) = powergenerate_monhour_CN-bbbb;
    powergenerate_monhour_reg(:,reg) = powergenerat_monhour(:,i) + powergenerate_monhour_reg(:,reg);% TWh/h
    powergenerate_monhour_CN2_2060_trans_storage1=powergenerate_monhour_CN2_2060_trans_storage;
%   power use efficiency in 2060
    for zz = 1:288
        if powergenerate_monhour_reg2_2060_trans_storage(zz,reg) <= powerdemand_monhour_reg(zz,reg) && powergenerate_monhour_CN2_2060_trans_storage1(zz,1) <= powerdemand_monhour_CN2060(zz,1)  && sum(sum(powergenerate_monhour_CN2_2060_trans_storage1))<=sum(sum(powerdemand_monhour_CN2060))
            ttt=1;
            % Intra-region consumption
            powergenerate_monhour_reg2_2060(zz,reg) = powergenerate_monhour_reg(zz,reg)*0.967;
            powergenerate_monhour_reg2_2060(zz,reg) = min(powerdemand_monhour_reg(zz,reg),powergenerate_monhour_reg2_2060(zz,reg));
            a1=-aaaa_reg(zz,reg) + powergenerate_monhour_reg2_2060(zz,reg);
            a1(a1<0)=0;
            powergenerate_monhour_CN2_2060(zz,1) = aaaa(zz,1)+a1;
             % considering transportation
            powergenerate_monhour_reg2_2060_trans(zz,reg) = powergenerate_monhour_reg(zz,reg)*0.967; 
            a1=-aaaa_reg_trans(zz,reg) + powergenerate_monhour_reg2_2060_trans(zz,reg);
            a1(a1<0)=0;
            powergenerate_monhour_CN2_2060_trans(zz,1) = aaaa_trans(zz,1)+a1;
            powergenerate_monhour_CN2_2060_trans(zz,1) = min(powerdemand_monhour_CN2060(zz,1),powergenerate_monhour_CN2_2060_trans(zz,1));
             % considering transportation and storage
            powergenerate_monhour_reg2_2060_trans_storage(zz,reg) = powergenerate_monhour_reg(zz,reg)*0.967; 
            a1=-aaaa_reg_trans_storage(zz,reg) + powergenerate_monhour_reg2_2060_trans_storage(zz,reg);
            a1(a1<0)=0;
            powergenerate_monhour_CN2_2060_trans_storage(zz,1) = aaaa_trans_storage(zz,1)+a1;
        end 
        if powergenerate_monhour_reg2_2060_trans_storage(zz,reg) > powerdemand_monhour_reg(zz,reg) && powergenerate_monhour_CN2_2060_trans_storage1(zz,1) <= powerdemand_monhour_CN2060(zz,1)  && sum(sum(powergenerate_monhour_CN2_2060_trans_storage1))<=sum(sum(powerdemand_monhour_CN2060))
            ttt=2;
            % Intra-region consumption
            powergenerate_monhour_reg2_2060(zz,reg) = powerdemand_monhour_reg(zz,reg);
            powergenerate_monhour_CN2_2060(zz,1) = aaaa(zz,1);
             % considering transportation
            powergenerate_monhour_reg2_2060_trans(zz,reg) = powerdemand_monhour_reg(zz,reg)+(powergenerate_monhour_reg(zz,reg)*0.967-powerdemand_monhour_reg(zz,reg))*0.983; 
            a1= - aaaa_reg_trans(zz,reg) + powergenerate_monhour_reg2_2060_trans(zz,reg);
            a1(a1<0)=0;
            powergenerate_monhour_CN2_2060_trans(zz,1) = aaaa_trans(zz,1)+a1;
            powergenerate_monhour_CN2_2060_trans(zz,1)=min(powergenerate_monhour_CN2_2060_trans(zz,1),powerdemand_monhour_CN2060(zz,1));
             % considering transportation and storage
            powergenerate_monhour_reg2_2060_trans_storage(zz,reg) = powerdemand_monhour_reg(zz,reg)+(powergenerate_monhour_reg(zz,reg)*0.967-powerdemand_monhour_reg(zz,reg))*0.983; 
            a1= -aaaa_reg_trans_storage(zz,reg) + powergenerate_monhour_reg2_2060_trans_storage(zz,reg);
            a1(a1<0)=0;
            powergenerate_monhour_CN2_2060_trans_storage(zz,1) = aaaa_trans_storage(zz,1)+a1;
        end
        if powergenerate_monhour_CN2_2060_trans_storage1(zz,1) > powerdemand_monhour_CN2060(zz,1) && sum(sum(powergenerate_monhour_CN2_2060_trans_storage1))<=sum(sum(powerdemand_monhour_CN2060))
           ttt=3;
            % Intra-region consumption
            powergenerate_monhour_reg2_2060(zz,reg) = aaaa_reg(zz,reg);
            powergenerate_monhour_CN2_2060(zz,1) = aaaa(zz,1);
             % considering transportation
            powergenerate_monhour_CN2_2060_trans(zz,1) = powerdemand_monhour_CN2060(zz,1);
             % considering transportation and storage
            powergenerate_monhour_CN2_2060_trans_storage(zz,1) =  aaaa_trans_storage(zz,1)+0.99*0.85*0.983*0.967*(powergenerat_monhour(zz,i));
            storage(zz,i) = storage_aaa(zz,i-1)+0.99*0.85*0.983*0.967*(powergenerat_monhour(zz,i));
        end
        if sum(sum(powergenerate_monhour_CN2_2060_trans_storage1))>sum(sum(powerdemand_monhour_CN2060))
            ttt=4;
            % Intra-region consumption
            powergenerate_monhour_reg2_2060(zz,reg) = aaaa_reg(zz,reg);
            powergenerate_monhour_CN2_2060(zz,1) = aaaa(zz,1);
             % considering transportation
            powerg = zeros(288,1);
            powergenerate_monhour_CN2_2060_trans(zz,1) = powerdemand_monhour_CN2060(zz,1);
             % considering transportation and storage
            powergenerate_monhour_CN2_2060_trans_storage(zz,1) = aaaa_trans_storage(zz,1)+0.99*0.85*0.983*0.967*(powergenerat_monhour(zz,i));
            storage(zz,i) = storage_aaa(zz,i-1)+0.99*0.85*0.983*0.967*powergenerat_monhour(zz,i);
        end
            powergenerate_monhour_CN2_20602(zz,i) = powergenerate_monhour_CN2_2060(zz,1)-aaaa(zz,1);
            powergenerate_monhour_CN2_2060_trans2(zz,i) = powergenerate_monhour_CN2_2060_trans(zz,1)-aaaa_trans(zz,1);
            powergenerate_monhour_CN2_2060_trans_storage2(zz,i) = powergenerate_monhour_CN2_2060_trans_storage(zz,1)-aaaa_trans_storage(zz,1);
    end
    fenzi(:,i)=powergenerate_monhour_CN2_2060-aaaa;

    fenzi_trans(:,i)=powergenerate_monhour_CN2_2060_trans-aaaa_trans;
    fenmu(:,i)=powergenerate_monhour_CN-bbbb;

    fenzi_trans_storage(:,i)=powergenerate_monhour_CN2_2060_trans_storage-aaaa_trans_storage;
    
    power_effi(:,i) = powergenerate_monhour_CN2_2060-aaaa;
    aaaa = powergenerate_monhour_CN2_2060;
    aaaa_reg = powergenerate_monhour_reg2_2060;
    aaaa_trans = powergenerate_monhour_CN2_2060_trans;
    aaaa_reg_trans = powergenerate_monhour_reg2_2060_trans;
    aaaa_trans_storage = powergenerate_monhour_CN2_2060_trans_storage;
    aaaa_reg_trans_storage = powergenerate_monhour_reg2_2060_trans_storage;
    
    bbbb = powergenerate_monhour_CN;
    powergenerate_aaa = powergenerate_monhour_CN2_2060_trans_storage;
    storage_aaa = storage;

    utilize_ratio2060(i,1) = sum(powergenerate_monhour_CN2_2060)./sum(powergenerate_monhour_CN);
    utilize_ratio2060_plant_alone(i,1) = (sum(sum(powergenerate_monhour_CN2_2060))-aa)./(sum(sum(powergenerate_monhour_CN))-bb);
    P(i,1) = (sum(sum(powergenerate_monhour_CN2_2060))-aa);
    P(i,2) = (sum(sum(powergenerate_monhour_CN))-bb);
    aa = sum(sum(powergenerate_monhour_CN2_2060));
    
    utilize_ratio2060_trans(i,1) = sum(powergenerate_monhour_CN2_2060_trans)./sum(powergenerate_monhour_CN);
    utilize_ratio2060_trans_plant_alone(i,1) = (sum(sum(powergenerate_monhour_CN2_2060_trans))-aa_trans)./(sum(sum(powergenerate_monhour_CN))-bb);
    P_trans(i,1) = (sum(sum(powergenerate_monhour_CN2_2060_trans))-aa_trans);
    P_trans(i,2) = (sum(sum(powergenerate_monhour_CN))-bb);
    aa_trans = sum(sum(powergenerate_monhour_CN2_2060_trans));
    
    utilize_ratio2060_trans_storage(i,1) = sum(powergenerate_monhour_CN2_2060_trans_storage)./sum(powergenerate_monhour_CN);
    utilize_ratio2060_trans_storage_plant_alone(i,1) = (sum(sum(powergenerate_monhour_CN2_2060_trans_storage))-aa_trans_storage)./(sum(sum(powergenerate_monhour_CN))-bb);
    P_trans_storage(i,1) = (sum(sum(powergenerate_monhour_CN2_2060_trans_storage))-aa_trans_storage);
    P_trans_storage(i,2) = (sum(sum(powergenerate_monhour_CN))-bb);
    aa_trans_storage = sum(sum(powergenerate_monhour_CN2_2060_trans_storage));
    bb = sum(sum(powergenerate_monhour_CN));
    i
end
utilize_ratio2060_trans_plant_alone(579)=utilize_ratio2060_trans_plant_alone(578)/2+utilize_ratio2060_trans_plant_alone(580)/2;
utilize_ratio2060_trans_storage_plant_alone(579)=utilize_ratio2060_trans_storage_plant_alone(578)/2+utilize_ratio2060_trans_storage_plant_alone(580)/2;

utilize_ratio2060_trans_plant_alone=max(utilize_ratio2060_plant_alone,utilize_ratio2060_trans_plant_alone);
storage1 = zeros(288,3844);
storage1(:,1)=storage(:,1);
for i = 2:3844
    if sum(storage(:,i))~=0
        storage1(:,i)= storage(:,i)-storage(:,i-1);
    end
end
storage_max_plant = max(storage1)'; % TWh/h
storage_year_plant = sum(storage1)'/288*8760; % TWh/year

save('G:\China C neutrality\Data\ANS\utilize_ratio2060_plant_alone.mat','utilize_ratio2060_plant_alone') 
save('G:\China C neutrality\Data\ANS\utilize_ratio2060_trans_plant_alone.mat','utilize_ratio2060_trans_plant_alone') 
save('G:\China C neutrality\Data\ANS\utilize_ratio2060_trans_storage_plant_alone.mat','utilize_ratio2060_trans_storage_plant_alone') 
save('G:\China C neutrality\Data\ANS\storage_max_plant.mat','storage_max_plant')   % TWh/h  
save('G:\China C neutrality\Data\ANS\storage_year_plant.mat','storage_year_plant')  % TWh/year

