tic
clear;
load('H:\China C neutrality\Data\CPfix_AC.mat');  % lines
% AC 1000 kV: 1:transmission distance (km); 2:transmission capacity (MW)	

rmb2us=1/6.8967; % RMB to USD2019
% fossilfuel_emissionfactor=0.85; % tCO2/MWh Brander, M. Electricity-specific Emission Factors for Grid Electricity 2011
fossilfuel_emissionfactor=0.783; 
CO2_C=0.2727;
% carbonprice=10; % USD2019/tonCO2
lifetime_power=25;
discount=0.05; % per year
discount1yr=0;
for t=1:lifetime_power
    discount1yr=discount1yr+1/(1+discount)^(t-1);
end
% PV:
OMratio_majorline_PV=0.01;
OMratio_substation_PV=0.01;
% wind
OMratio_majorline_wind=0.03;
OMratio_substation_wind=0.03;
load('H:\China C neutrality\PV_power potential\ANS_PV1\tranmission_lines.dat','-mat');  % lines
% Column 15: Distance between two stations
% Column 16: Type of UHV, 1 denotes 1000kV alternating current power
% station,2 denotes ±1100kV direct  current power station, 3 denotes ±800kV direct  current power station
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
lines(numlines+1:end,:)=[];

for i = 1:numlines
    type(i,1) = lines(i,16);
    if type(i,1) == 1
        CP_fix(i,1) = interp1(CPfix_AC(:,1), CPfix_AC(:,2),lines(i,15),'linear'); % MW
        cost_line(i,1) = 670785; % $/km
        cost_sta(i,1) = 41*10^3; % $ MW-1
    end
    if type(i,1) == 2
        CP_fix(i,1) = 12000; %MW
        cost_line(i,1) = 800383; % $/km
        cost_sta(i,1) = 92*10^3; % $ MW-1
    end
    if type(i,1) == 3
        CP_fix(i,1) = 8000; %MW
        cost_line(i,1) = 732220; % $/km
        cost_sta(i,1) = 82*10^3; % $ MW-1
    end
end

load('H:\China C neutrality\PV_power potential\ANS_PV1\optpowerunit_PV.mat'); % 
load('H:\China C neutrality\PV_power potential\ANS_PV1\powerunit_IX_PV.mat'); % 
load('H:\China C neutrality\PV_power potential\ANS_PV1\tranmission_lines_IX.mat');  % lines_IX
lines_IX(size(optpowerunit_PV,1)+1:end,:)=[];
lines_IX_PV = lines_IX;
load('H:\China C neutrality\PV_power potential\ANS_PV1\powerunit_num_IX_PV.mat');  % 
optpowerunit_PV(:,35) = 1;
optpowerunit_PV(:,40) = powerunit_IX_PV; % Power plant ID

load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\optpowerunit_onshorewind.mat'); % 
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit_IX_onshorewind.mat'); % 
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\tranmission_lines_IX.mat');  % lines_IX
lines_IX(size(optpowerunit_onshorewind,1)+1:end,:)=[];
lines_IX_onshorewind = lines_IX;
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit_num_IX_onshorewind.mat');  % 
optpowerunit_onshorewind(:,35) = 2;
optpowerunit_onshorewind(:,40) = powerunit_IX_onshorewind; % Power plant ID

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
optpowerunit_offshorewind(:,40) = powerunit_IX_offshorewind; % Power plant ID

optpowerunit = [optpowerunit_PV;optpowerunit_onshorewind;optpowerunit_offshorewind];
lines_IX = [lines_IX_PV(:,1:15);lines_IX_onshorewind;lines_IX_offshorewind];
powerunit_num_IX = [powerunit_num_IX_PV;powerunit_num_IX_onshorewind;powerunit_num_IX_offshorewind];
[B,IX]=sort(optpowerunit(:,20),1);
numpowerunit = size(optpowerunit,1);
for i=1:numpowerunit
    i2=IX(i);
    powerunit_IX(i,1)=i2;
    optpowerunit_IX(i,1:40)=optpowerunit(i2,1:40); % lat lon
    lines_IX_IX(i,1:15)=lines_IX(i2,1:15); % lat lon
    powerunit_num_IX_IX(i,1:4)=powerunit_num_IX(i2,1:4); % lat lon
end

line_IX_all = [lines(:,1:15);lines_IX_IX];
[m,n] = find(line_IX_all(:,9)==0);
line_IX_all(m,9) = 1362; %offshorewind, Shenzhen
%%
rate_energy=0.02; % per year
load('H:\China C neutrality\Data\powerdemand_pro2060_ele_0424.mat')  % power demand in 2060 with 58% and 55.6% electrification，substracting bioenergy, nuclear, hydro, hydrogen energy, TWh/year

load('H:\China C neutrality\Data\ID_Pro_CN.mat');  % ID_Pro_CN(4800x1950)
powerdemand_pro2060_34 = zeros(34,1);
for i = 1:max(max(ID_Pro_CN))
        [m,n] = find(ID_Pro_CN==i);
        powerdemand_pro2060_34(i,1)=sum(sum(powerdemand_pro2060_2(sub2ind(size(powerdemand_pro2060_2), m, n))));
        clear m
        clear n
        i
end
powerdemand_pro2060_34(powerdemand_pro2060_34<0)=0;

%%
E_120_2060 = powerdemand_pro2060_2*1000/365; % TWh/year-> GWh/grid/day
E_120 = E_120_2060./((1+rate_energy)^(2060-2019)); % GWh/grid/day 2019年

load('H:\China C neutrality\Data\ID_Pro_CN.mat');  % ID_Pro_CN(4800x1950)
load('H:\China C neutrality\Data\ID_County_CN.mat'); % county_CN(4800x1950)
land_sink=load('H:\China C neutrality\Data\land_sink_sr2.txt'); % land_sink(72x46)
% LC_CN(72x46) gC/m2/yr land carbon sink from atmospheric inversion
% Wang et al.Large Chinese land carbon sink estimated from atmospheric carbon dioxide data. Nature, 2021

Rearth    =  6371.3;      % km average radium of the earth
E_120cy=zeros(2373,2); % electricity demand GWh/day from 2019 to 2060
e0=zeros(8,5); % Mton C 1 land carbon sink; 2 carbon
gridarea=zeros(4800,1950);
regions = [2,3,2,5,6,6,7,6,3,1,6,4,1,1,4,2,2,4,0,3,5,5,7,5,2,2,3,3,6,5,7,7,2,7];
for i=1:4800
    for j=1:1950
        gridarea(i,j)=abs(Rearth^2*(sin(((55-(i-1)/120+1/240)+1/120)*pi/180)-sin((55-(i-1)/120+1/240)*pi/180))*1/30*pi/180); %km2
        prv=ID_Pro_CN(i,j);
        if prv>=1 && prv<=34
            reg=regions(prv);
            if reg>=1 && reg<=7
                e0(reg,1)=e0(reg,1)+E_120(i,j); % GWh/day
                e0(reg,5)=e0(reg,5)+land_sink(floor((73+(j-1)/30+180)/5)+1,floor((55-(i-1)/120+90)/4)+1)*gridarea(i,j)*1e-6; % land carbon sink Mton C/year
            end
        end
        cy=county_CN(i,j);
        if cy>=1 && cy<=2373
            E_120cy(cy,1)=E_120cy(cy,1)+E_120(i,j); % GWh/day
        end
    end
end
rate_energy=0.02; % per year
e0(8,:)=sum(e0(1:7,:),1); % 1 electricity demand; 5 land carbon sink (negative value) Mton C/year
e0(:,1)=e0(:,1).*(1+rate_energy)^(2060-2019)*0.365; % electricity demand by region in 2060 	TWh/year
E_120cy(:,2)=E_120cy(:,1).*(1+rate_energy)^(2060-2019); % electricity demand by county in 2060
emission2060=e0(8,1)*fossilfuel_emissionfactor*CO2_C+e0(8,5); % CO2 emission from fossil fuel offset by land sink Mt C / year

e0_PV = e0;
e0_onshorewind = e0;
e0_offshorewind = e0;

[B,IX]=sort(optpowerunit_IX(:,20),1);
linescapacity=zeros(numlines,1); % capacity of electricity transmissions
linescapacity0=linescapacity;
linescapacity_plant=zeros(numlines,3844);
etrans=zeros(7,7); % nenux of electricity transmissions: etrans(r1,r2) is from r1 to r2
etrans_PV=zeros(7,7); %
etrans_onshorewind=zeros(7,7); %
etrans_offshorewind=zeros(7,7); %
topon=[1 2 2 4 4];
load('H:\China C neutrality\ANS\utilize_ratio2060_trans_storage_plant_alone.mat')  % 288*1 考虑用电效率、区域间运输和储存
optpowerunit_IX(:,1)=optpowerunit_IX(:,1).*utilize_ratio2060_trans_storage_plant_alone;
for i=1:numpowerunit
    i2=IX(i);
%     powerunit_IX(i,1)=i2;
    % electricity transmission
    cy=line_IX_all(numlines+i2,9); % county of power unit
    reg=line_IX_all(numlines+i2,11); % region of power unit 1-7
    e0(reg,2)=e0(reg,2)+optpowerunit_IX(i2,1); % electricity produced by the region TWh / year
    if optpowerunit_IX(i2,35)==1
        e0_PV(reg,2)=e0_PV(reg,2)+optpowerunit_IX(i2,1); % electricity produced by the region TWh / year
    else if optpowerunit_IX(i2,35)==2
            e0_onshorewind(reg,2)=e0_onshorewind(reg,2)+optpowerunit_IX(i2,1); % electricity produced by the region TWh / year
        else if optpowerunit_IX(i2,35)==3
            e0_offshorewind(reg,2)=e0_offshorewind(reg,2)+optpowerunit_IX(i2,1); % electricity produced by the region TWh / year
            end
        end
    end
    powercy=E_120cy(cy,2)/topon(powerunit_num_IX_IX(i2,4))*0.365; % electricity demand in this county GWh/day -> TWh/year    
    if (e0(reg,4)+powercy)>=e0(reg,1)
        powercy=0; % the electricity demand has been 100% met by transmission
    end
    if optpowerunit_IX(i2,1)<=powercy
        e0(reg,3)=e0(reg,3)+optpowerunit_IX(i2,1); % electricity used by the county TWh / year
        e0(reg,4)=e0(reg,4)+optpowerunit_IX(i2,1); % electricity used by the region TWh / year
        etrans(reg,reg)=etrans(reg,reg)+optpowerunit_IX(i2,1); % electricity from reg to reg TWh / year
        if optpowerunit_IX(i2,35)==1
            e0_PV(reg,3)=e0_PV(reg,3)+optpowerunit_IX(i2,1); % electricity used by the county TWh / year
            e0_PV(reg,4)=e0_PV(reg,4)+optpowerunit_IX(i2,1); % electricity used by the region TWh / year
            etrans_PV(reg,reg)=etrans_PV(reg,reg)+optpowerunit_IX(i2,1); % electricity from reg to reg TWh / year
        else if optpowerunit_IX(i2,35)==2
                e0_onshorewind(reg,3)=e0_onshorewind(reg,3)+optpowerunit_IX(i2,1); % electricity used by the county TWh / year
                e0_onshorewind(reg,4)=e0_onshorewind(reg,4)+optpowerunit_IX(i2,1); % electricity used by the region TWh / year
                etrans_onshorewind(reg,reg)=etrans_onshorewind(reg,reg)+optpowerunit_IX(i2,1); % electricity from reg to reg TWh / year
            else if optpowerunit_IX(i2,35)==3
                    e0_offshorewind(reg,3)=e0_offshorewind(reg,3)+optpowerunit_IX(i2,1); % electricity used by the county TWh / year
                    e0_offshorewind(reg,4)=e0_offshorewind(reg,4)+optpowerunit_IX(i2,1); % electricity used by the region TWh / year
                    etrans_offshorewind(reg,reg)=etrans_offshorewind(reg,reg)+optpowerunit_IX(i2,1); % electricity from reg to reg TWh / year
                end
            end
        end
    else
        e0(reg,3)=e0(reg,3)+powercy; % electricity used by the county TWh / year
        e0(reg,4)=e0(reg,4)+powercy; % electricity used by the region TWh / year
        etrans(reg,reg)=etrans(reg,reg)+powercy; % electricity from reg to reg TWh / year
        if optpowerunit_IX(i2,35)==1
            e0_PV(reg,3)=e0_PV(reg,3)+powercy; % electricity used by the county TWh / year
            e0_PV(reg,4)=e0_PV(reg,4)+powercy; % electricity used by the region TWh / year
            etrans_PV(reg,reg)=etrans_PV(reg,reg)+powercy; % electricity from reg to reg TWh / year
        else if optpowerunit_IX(i2,35)==2
                e0_onshorewind(reg,3)=e0_onshorewind(reg,3)+powercy; % electricity used by the county TWh / year
                e0_onshorewind(reg,4)=e0_onshorewind(reg,4)+powercy; % electricity used by the region TWh / year
                etrans_onshorewind(reg,reg)=etrans_onshorewind(reg,reg)+powercy; % electricity from reg to reg TWh / year
            else if optpowerunit_IX(i2,35)==3
                    e0_offshorewind(reg,3)=e0_offshorewind(reg,3)+powercy; % electricity used by the county TWh / year
                    e0_offshorewind(reg,4)=e0_offshorewind(reg,4)+powercy; % electricity used by the region TWh / year
                    etrans_offshorewind(reg,reg)=etrans_offshorewind(reg,reg)+powercy; % electricity from reg to reg TWh / year
                end
            end
        end
        % transmission of energy by the first major line
        powertrans=optpowerunit_IX(i2,1)-powercy; %
        lineid=line_IX_all(numlines+i2,6); % id of line endpoint
        linescapacity_plant(lineid,i2)=linescapacity_plant(lineid,i2)+optpowerunit_IX(i2,30)*powertrans/optpowerunit_IX(i2,1); % transmission capacity MW
        linescapacity(lineid,1)=linescapacity(lineid,1)+optpowerunit_IX(i2,30)*powertrans/optpowerunit_IX(i2,1); % transmission capacity MW
        reg2=line_IX_all(lineid,14); % region of endpoint 1-7
        if e0(reg2,1)>=(e0(reg2,4)+powertrans)
            etrans(reg,reg2)=etrans(reg,reg2)+powertrans; % electricity from reg to reg2 TWh / year
            e0(reg2,4)=e0(reg2,4)+powertrans; % electricity used by the region TWh / year
            if optpowerunit_IX(i2,35)==1
                etrans_PV(reg,reg2)=etrans_PV(reg,reg2)+powertrans; % electricity from reg to reg2 TWh / year
                e0_PV(reg2,4)=e0_PV(reg2,4)+powertrans; % electricity used by the region TWh / year
            else if optpowerunit_IX(i2,35)==2
                    etrans_onshorewind(reg,reg2)=etrans_onshorewind(reg,reg2)+powertrans; % electricity from reg to reg2 TWh / year
                    e0_onshorewind(reg2,4)=e0_onshorewind(reg2,4)+powertrans; % electricity used by the region TWh / year
                else if optpowerunit_IX(i2,35)==3
                        etrans_offshorewind(reg,reg2)=etrans_offshorewind(reg,reg2)+powertrans; % electricity from reg to reg2 TWh / year
                        e0_offshorewind(reg2,4)=e0_offshorewind(reg2,4)+powertrans; % electricity used by the region TWh / year
                    end
                end
            end
        else
            % electricity exceedes the demand
            etrans(reg,reg2)=etrans(reg,reg2)+e0(reg2,1)-e0(reg2,4); % electricity from reg to reg2 TWh / year
            powertrans2=(e0(reg2,4)+powertrans)-e0(reg2,1); % electricity exceeding demand TWh / year
            aaa = e0(reg2,1)-e0(reg2,4);
            e0(reg2,4)=e0(reg2,1);
            if optpowerunit_IX(i2,35)==1
                etrans_PV(reg,reg2)=etrans_PV(reg,reg2)+aaa; % electricity from reg to reg2 TWh / year
                e0_PV(reg2,4)=e0_PV(reg2,4)+aaa;
            else if optpowerunit_IX(i2,35)==2
                    etrans_onshorewind(reg,reg2)=etrans_onshorewind(reg,reg2)+aaa; % electricity from reg to reg2 TWh / year
                    e0_onshorewind(reg2,4)=e0_onshorewind(reg2,4)+aaa;
                else if optpowerunit_IX(i2,35)==3
                        etrans_offshorewind(reg,reg2)=etrans_offshorewind(reg,reg2)+aaa; % electricity from reg to reg2 TWh / year
                        e0_offshorewind(reg2,4)=e0_offshorewind(reg2,4)+aaa;
                    end
                end
            end
           
            % find transmission of electricity by the same major line
            idx=find(line_IX_all(1:numlines,5)==line_IX_all(lineid,5));
            reg3=line_IX_all(idx(end),14); % region of endpoint 1-7
            if reg3~=reg2 && (e0(reg3,4)+powertrans2)<=e0(reg3,1)
                % transmission of electricity to the final destination of line
                etrans(reg,reg3)=etrans(reg,reg3)+powertrans2; % electricity from reg to reg3 TWh / year
                e0(reg3,4)=e0(reg3,4)+powertrans2; % electricity used by reg 3 TWh / year
                if optpowerunit_IX(i2,35)==1
                    etrans_PV(reg,reg3)=etrans_PV(reg,reg3)+powertrans2; % electricity from reg to reg3 TWh / year
                    e0_PV(reg3,4)=e0_PV(reg3,4)+powertrans2; % electricity used by reg 3 TWh / year
                else if optpowerunit_IX(i2,35)==2
                        etrans_onshorewind(reg,reg3)=etrans_onshorewind(reg,reg3)+powertrans2; % electricity from reg to reg3 TWh / year
                        e0_onshorewind(reg3,4)=e0_onshorewind(reg3,4)+powertrans2; % electricity used by reg 3 TWh / year
                    else if optpowerunit_IX(i2,35)==3
                            etrans_offshorewind(reg,reg3)=etrans_offshorewind(reg,reg3)+powertrans2; % electricity from reg to reg3 TWh / year
                            e0_offshorewind(reg3,4)=e0_offshorewind(reg3,4)+powertrans2; % electricity used by reg 3 TWh / year
                        end
                    end
                end
                
                linescapacity(lineid:idx(end),1)=linescapacity(lineid:idx(end),1)+optpowerunit_IX(i2,30)*powertrans2/optpowerunit_IX(i2,1); % transmission capacity MW
                linescapacity_plant(lineid:idx(end),i2)=linescapacity_plant(lineid:idx(end),i2)+optpowerunit_IX(i2,30)*powertrans2/optpowerunit_IX(i2,1); % transmission capacity MW
            else
                % transmission of electricity from reg2 to the closest and largest center reg4
                [B2,IX2]=sort(e0(1:7,4)-e0(1:7,1),1);
                reg4=IX2(1);
                for j=1:numlines
                    if line_IX_all(j,14)==reg4
                        j2=j;
                        if line_IX_all(j,11)==reg2
                            j2=j;
                        end
                    end
                end
                
                etrans(reg,reg4)=etrans(reg,reg4)+powertrans2; % electricity from reg to reg4 TWh / year
                e0(reg4,4)=e0(reg4,4)+powertrans2; % electricity used by reg4 TWh / year
                if optpowerunit_IX(i2,35)==1
                    etrans_PV(reg,reg4)=etrans_PV(reg,reg4)+powertrans2; % electricity from reg to reg3 TWh / year
                    e0_PV(reg4,4)=e0_PV(reg4,4)+powertrans2; % electricity used by reg 3 TWh / year
                else if optpowerunit_IX(i2,35)==2
                        etrans_onshorewind(reg,reg4)=etrans_onshorewind(reg,reg4)+powertrans2; % electricity from reg to reg3 TWh / year
                        e0_onshorewind(reg4,4)=e0_onshorewind(reg4,4)+powertrans2; % electricity used by reg 3 TWh / year
                    else if optpowerunit_IX(i2,35)==3
                            etrans_offshorewind(reg,reg4)=etrans_offshorewind(reg,reg4)+powertrans2; % electricity from reg to reg3 TWh / year
                            e0_offshorewind(reg4,4)=e0_offshorewind(reg4,4)+powertrans2; % electricity used by reg 3 TWh / year
                        end
                    end
                end
                linescapacity_plant(j2,i2)=linescapacity_plant(j2,i2)+optpowerunit_IX(i2,30)*powertrans2/optpowerunit_IX(i2,1); % transmission capacity MW
                linescapacity(j2,1) = linescapacity(j2,1)+optpowerunit_IX(i2,30)*powertrans2/optpowerunit_IX(i2,1);  % transmission capacity MW
            end
        end
    end
    for j=1:7
        e0(j,4)=min(e0(j,1),e0(j,4));
    end
    e0(8,1:end)=sum(e0(1:7,1:end),1);
    if optpowerunit_IX(i2,35)==1
        for j=1:7
            e0_PV(j,4)=min(e0_PV(j,1),e0_PV(j,4));
        end
        e0_PV(8,1:end)=sum(e0_PV(1:7,1:end),1);
    else if optpowerunit_IX(i2,35)==2
            for j=1:7
                e0_onshorewind(j,4)=min(e0_onshorewind(j,1),e0_onshorewind(j,4));
            end
            e0_onshorewind(8,1:end)=sum(e0_onshorewind(1:7,1:end),1);
        else if optpowerunit_IX(i2,35)==3
                for j=1:7
                    e0_offshorewind(j,4)=min(e0_offshorewind(j,1),e0_offshorewind(j,4));
                end
                e0_offshorewind(8,1:end)=sum(e0_offshorewind(1:7,1:end),1);
            end
        end
    end
     i
    eee_reg_plant(:,i) = e0(:,4); % electricity used by reg 3 TWh / year
end

% fixed capacity of 8000 MW, 12000 MW for ±800 kV DC and ±1100kV DC
% $732220, $800383 and $670785 km-1 for ±800 kV DC, ±1100kV DC and 1000 kV AC, respectively
% $82, $92 and $41 kW-1 for ±800 kV DC, ±1100kV DC and 1000 kV AC, respectively
% CP_fix %MW
% cost_line % $/km
% cost_sta % $ MW-1
% linescapacity %MW

% Column 1-15 are 1000kV alternating current power stations, Column 16 is ±1100kV direct  current power station,  Column 17-39 are ±800kV direct  current power stations 
num_line_0 = lines(:,16);
num_line_0(num_line_0==1) = 0.1;
num_line_0(num_line_0==2) = 1;
num_line_0(num_line_0==3) = 1;
num_line_0(num_line_0==0.1) = 2;

cost_l_all = (linescapacity.*cost_sta+lines(:,15).*cost_line.*max(ceil(linescapacity./CP_fix)-num_line_0,0))/10^6; % million $

for i = 1:numlines
    cost_l_plant(i,:) = cost_l_all(i,1)*(linescapacity_plant(i,:)./sum(linescapacity_plant(i,:))); % million $
end
cost_l_plant(find(isnan(cost_l_plant)==1))=0;
cost_trans_IX = sum(cost_l_plant)';% million $
CP_trans_IX = sum(linescapacity_plant); % MW

save('H:\China C neutrality\ANS\cost_trans_IX2.mat','cost_trans_IX')% million $
save('H:\China C neutrality\ANS\CP_trans_IX2.mat','CP_trans_IX')% million $



