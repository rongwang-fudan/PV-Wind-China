tic
clear;
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\tranmission_lines.dat','-mat'); % lines
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\ID_dist_offshore120.mat'); % $/W
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\S_offshorewind.dat','-mat'); 
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\costunits_offshorewind.dat','-mat');
% 1 pro; 2 Powerunit_ID; 3 CP,kW; 4 Cost,million dollar; 5 kwh/year; 
% 6 S km2; 7 mean wind speeed,m/s 
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\pro_offshore120.mat'); % pro ID
% 3:Fujian; 5:Guangdong; 6: Guangxi; 8:Hainan; 9:Hebei; 17:Jiangsu; 
% 18:Liaoning; 25: Shandong; 26:Shanghai; 28: Tianjin; 32: Yunnan; 33:Zhejiang; 
load('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\ID_dist_offshore120.mat'); % $/W
SRTM30_CN=imread('G:\China C neutrality\Data\SRTM30_CN120.tif');
SRTM30_CN = double(SRTM30_CN);
SRTM30_CN = SRTM30_CN(1:4800,:);
module_price = 1.008; % USD/W
rmb2us=1/6.8967; % RMB to USD2019
lifetime_power=40;
discount=0.07; % per year
discount1yr=0;
for t=1:lifetime_power
    discount1yr=discount1yr+1/(1+discount)^(t-1);
end
degration = 0.0; 
degrat1yr=0;
for t=1:lifetime_power
    degrat1yr=degrat1yr+(1-degration)^(t-1)/(1+discount)^(t-1);
end
OMratio_substation=1/lifetime_power*1.2;
fossilfuel_emissionfactor=0.827; % tCO2/MWh Brander
CO2_C=0.2727;
carbonprice=10; % USD2019/tonCO2
pro_ID = [3 5 6 8 9 17 18 25 26 28 33];
offshorewind_power_lcoemin = zeros(4800,7800);
unitid_lcoe = zeros(4800,7800);
numpowerunit=11;
for i = 1:numpowerunit
    lcoemin=1000;
    [m,n] = find(costunits_offshorewind(:,1)==pro_ID(i));
    lcoeunit=zeros(size(m,1),8);
    wind_speed=zeros(size(m,1),1);
    for j=m(1):m(end)
        j1=j-m(1)+1;
        if j-m(1)>0
            lcoeunit(j1,1)=lcoeunit(j1-1,1)+costunits_offshorewind(j,5)/10^9; % electricity TWh / year
            lcoeunit(j1,2)=lcoeunit(j1-1,2)+costunits_offshorewind(j,4); % Cost,million dollar
            lcoeunit(j1,3)=lcoeunit(j1-1,3)+costunits_offshorewind(j,3); % CP,kW
            lcoeunit(j1,4)=lcoeunit(j1-1,4)+costunits_offshorewind(j,6); % S km2;
            lcoeunit(j1,5)=lcoeunit(j1-1,5)-costunits_offshorewind(j,5)/10^9*fossilfuel_emissionfactor*CO2_C; % abated annual CO2 emission from fossil fuel Mton C / year
            lcoeunit(j1,9)=lcoeunit(j1-1,9)+costunits_offshorewind(j,3)*1000*module_price/10^6; % module cost, million USD
            wind_speed(j1,1)=wind_speed(j1-1,1)+costunits_offshorewind(j,7); % mean wind speeed,m/s 
        else
            lcoeunit(j1,1)=costunits_offshorewind(j,5)/10^9; % electricity TWh / year
            lcoeunit(j1,2)=costunits_offshorewind(j,4); % Cost,million dollar
            lcoeunit(j1,3)=costunits_offshorewind(j,3); % CP,kW
            lcoeunit(j1,4)=costunits_offshorewind(j,6); % S km2;
            lcoeunit(j1,5)=-costunits_offshorewind(j,5)/10^9*fossilfuel_emissionfactor*CO2_C; % abated annual CO2 emission from fossil fuel Mton C / year
            lcoeunit(j1,9)=costunits_offshorewind(j,3)*1000*module_price/10^6; % module cost, million USD
            wind_speed(j1,1)=costunits_offshorewind(j,7); % mean wind speeed,m/s 
        end
        lcoeunit(j1,6)=lcoeunit(j1,2)*(1+OMratio_substation*discount1yr); % cost
        lcoeunit(j1,10)=lcoeunit(j1,9)*(1+OMratio_substation*discount1yr); % cost
        
        % LCoE
        lcoeunit(j1,8)=(lcoeunit(j1,6)+lcoeunit(j1,7))/lcoeunit(j1,1)/degrat1yr/1000; % LCoE million USD2019/TWh->USD2019/kWh        
        if lcoeunit(j1,8)<lcoemin
            jopt=j1;
            lcoemin=lcoeunit(j1,8);
        end
    end
        pro_ID_minlcoe(i,1)=jopt;
        jj=jopt;
        powers(i,1:10)=lcoeunit(jj,1:10); % 20 for LCoE USD2019/kWh
        wind_speed1(i,1) = wind_speed(jj,1)/jj;
        
        [m,n]=find(pro_offshore120==pro_ID(i) & ID_dist_offshore120<=pro_ID_minlcoe(i,1) & S_offshorewind~=0 & SRTM30_CN<=0);
        offshorewind_power_lcoemin(sub2ind(size(offshorewind_power_lcoemin), m, n))= sum(sum(powers(i,1)));%pro_offshore120(sub2ind(size(pro_offshore120), m, n));
        unitid_lcoe(sub2ind(size(unitid_lcoe), m, n))= pro_offshore120(sub2ind(size(pro_offshore120), m, n));%pro_offshore120(sub2ind(size(pro_offshore120), m, n));
end
lcoe_offshorewind=powers(:,8);
save('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\lcoe_offshorewind.mat','lcoe_offshorewind'); % USD2019/kWh
save('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\unitid_lcoe.mat','unitid_lcoe'); % USD2019/kWh
save('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\offshorewind_power_lcoemin.mat','offshorewind_power_lcoemin'); % power ID
%%
regions = [2,3,2,5,6,6,7,6,3,1,6,4,1,1,4,2,2,4,0,3,5,5,7,5,2,2,3,3,6,5,7,7,2,7];
CO2=zeros(8,1); % Mton C 1 land carbon sink; 2 carbon
load('G:\China C neutrality\Data\powerunit_offshorewind.mat'); % lat-y, lon-x, province
powerunit=powerunit_offshorewind;
powerunit(:,1)=floor((55-powerunit(:,1))*120)+1; % conversion of lat to y
powerunit(:,2)=floor((powerunit(:,2)-73)*30)+1; % conversion of lon to x
for i=1:11
    prv=pro_ID(i);
    if prv>=1 && prv<=34
        reg=regions(prv);
        CO2(reg,1)=CO2(reg,1)+powers(i,5); % abated annual CO2 emission from fossil fuel Mton C / year
    end
end
CO2(8,:)=sum(CO2(1:7,:),1); 

%%
% Energy Consumption 3
load('G:\China C neutrality\Data\powerdemand_pro2060_2.mat')  
rate_energy=0.02; % per year
E_120_2060 = powerdemand_pro2060_2*1000/365; % TWh/year-> GWh/grid/day
E_120_2 = E_120_2060./((1+rate_energy)^(2060-2019)); % GWh/grid/day 2019年
clear E_120
for i = 1:4800
    for j = 1:1950
        E_120(i,j*4-3:j*4) = E_120_2(i,j)/4;
    end
    i
end
clear E_120_2
load('G:\China C neutrality\Data\ID_Pro_CN.mat');  % ID_Pro_CN(4800x1950)
load('G:\China C neutrality\Data\ID_County_CN.mat'); % county_CN(4800x1950)
for i = 1:4800
    for j = 1:1950
        ID_Pro_CN1(i,j*4-3:j*4) = ID_Pro_CN(i,j);
        county_CN1(i,j*4-3:j*4) = county_CN(i,j);
    end
end
land_sink=load('G:\China C neutrality\Data\land_sink_sr2.txt'); % land_sink(72x46)
% LC_CN(72x46) gC/m2/yr land carbon sink from atmospheric inversion
% Wang et al.Large Chinese land carbon sink estimated from atmospheric carbon dioxide data. Nature, 2021

Rearth    =  6371.3;      % km average radium of the earth
E_120pro=zeros(34,2); % electricity demand GWh/day from 2019 to 2060
e0=zeros(8,5); % Mton C 1 land carbon sink; 2 carbon
gridarea=zeros(4800,1950);
regions = [2,3,2,5,6,6,7,6,3,1,6,4,1,1,4,2,2,4,0,3,5,5,7,5,2,2,3,3,6,5,7,7,2,7];
for i=1:4800
    for j=1:7800
        gridarea(i,j)=abs(Rearth^2*(sin(((55-(i-1)/120+1/240)+1/120)*pi/180)-sin((55-(i-1)/120+1/240)*pi/180))*1/30*pi/180); %km2
        prv=ID_Pro_CN1(i,j);
        if prv>=1 && prv<=34
            reg=regions(prv);
            if reg>=1 && reg<=7
                e0(reg,1)=e0(reg,1)+E_120(i,j); % GWh/day
            end
        end
        pro=ID_Pro_CN1(i,j);
        if pro>=1 && pro<=34
            E_120pro(pro,1)=E_120pro(pro,1)+E_120(i,j); % GWh/day
        end
    end
end
rate_energy=0.02; % per year
e0(8,:)=sum(e0(1:7,:),1); % 1 electricity demand; 5 land carbon sink (negative value) Mton C/year
e0(:,1)=e0(:,1).*(1+rate_energy)^(2060-2019)*0.365; % electricity demand by region in 2060 	TWh/year
E_120pro(:,2)=E_120pro(:,1).*(1+rate_energy)^(2060-2019); % electricity demand by province in 2060
emission2060=e0(8,1)*fossilfuel_emissionfactor*CO2_C+e0(8,5); % CO2 emission from fossil fuel offset by land sink Mt C / year

[B,IX]=sort(powers(:,8),1); 
optpowerunit=zeros(numpowerunit,36);
lines_IX=zeros(size(lines,1),size(lines,2));
powerunit_num_IX=zeros(size(powerunit,1),size(powerunit,2));
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
linescapacity=zeros(numlines,1); % capacity of electricity transmissions
linescapacity0=linescapacity;
etrans=zeros(7,7); % nenux of electricity transmissions: etrans(r1,r2) is from r1 to r2
topon=[1 2 2 4 4];
for i=1:numpowerunit
    i2=IX(i);
    powerunit_IX(i,1)=i2;
    powerunit_num_IX(i,:) = powerunit(i2,:);
    lines_IX(i,:)=lines(numlines+i2,:);
    optpowerunit(i,1:10)=powers(i2,1:10); % lat lon
    optpowerunit(i,31)=i2; % id in powerunit
    % national carbon budget in 2060
    optpowerunit(i,32)=optpowerunit(i,5);  % net CO2 emissions Mt C / year = (Abatement)
    if i>1
        optpowerunit(i,33)=optpowerunit(i-1,33)+optpowerunit(i,32); % net CO2 emissions Mt C / year
        optpowerunit(i,34)=optpowerunit(i-1,34)+optpowerunit(i,1); % electricity TWh / year
    else
        optpowerunit(i,33)=emission2060+optpowerunit(i,32); % net CO2 emissions Mt C / year = FF+(Abatement-Landsink)
        optpowerunit(i,34)=optpowerunit(i,1); % electricity TWh / year
    end
    % electricity transmission
    pro=lines(numlines+i2,10); % province of power unit
    reg=lines(numlines+i2,11); % region of power unit 1-7
end
optpowerunit_offshorewind = optpowerunit(:,1:34);
save('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\optpowerunit_offshorewind.mat','optpowerunit_offshorewind'); % power ID
powerunit_IX_offshorewind = powerunit_IX;
save('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\powerunit_IX_offshorewind.mat','powerunit_IX_offshorewind'); % power ID
save('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\tranmission_lines_IX.mat','lines_IX'); %
powerunit_num_IX_offshorewind = powerunit_num_IX;
save('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\powerunit_num_IX_offshorewind.mat','powerunit_num_IX_offshorewind'); %
