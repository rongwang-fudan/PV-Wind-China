%%
tic
clear;
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\wind_12to18_power.dat','-mat');
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit.dat','-mat');
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\unitid.dat','-mat');
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\unitid_lcoe.dat','-mat');
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\tranmission_lines.dat','-mat'); % lines
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\costunits.dat','-mat');
numpowerunit=size(powerunit,1);
rmb2us=1/6.8967; % RMB to USD2019
moduleprice=2.7*rmb2us; % USD2019/Wp  Module
consprice=0.85*rmb2us; % USD2019/Wp 
gridprice=0.17*rmb2us; % USD2019/Wp
otherprice=0.11*rmb2us; % USD2019/Wp
moduleprice2= moduleprice+consprice+gridprice+ otherprice; %  USD2019/Wp 
ratio_module = moduleprice/moduleprice2;

lifetime_power=25;
discount=0.05; % per year
discount1yr=0;
for t=1:lifetime_power
    discount1yr=discount1yr+1/(1+discount)^(t-1);
end
degration = 0.0;
degrat1yr=0;
for t=1:lifetime_power
    degrat1yr=degrat1yr+(1-degration)^(t-1)/(1+discount)^(t-1);
end

OMratio_majorline=0.03; %1/lifetime_power*1.2;
OMratio_substation=0.03; %1/lifetime_power*1.2;
fossilfuel_emissionfactor=0.783; 
CO2_C=0.2727;
carbonprice=0; % USD2019/tonCO2
powers=zeros(numpowerunit,30);

for i=1:numpowerunit
    idx=find(costs(:,1)==i);
    lcoeunit=zeros(size(idx,1),20);
    lcoeunit1=zeros(size(idx,1),1);
    SRR=zeros(size(idx,1),2);
    jopt=0;
    lcoemin=1000;
    lcoemin1=1000;
    for j=1:size(idx,1)
        electricity=costs(idx(j),3); % electricity TWh / year
        capacity=costs(idx(j),11); % capacity potential MW
        if j>1
            CPPP(j,1)=CPPP(j-1,1)+capacity; % capacity potential MW
            SRR(j,1)=SRR(j,1)+wind(idx(j),1);
            SRR(j,2)=SRR(j,2)+wind(idx(j),2);            
            lcoeunit(j,1)=lcoeunit(j-1,1)+electricity; % electricity TWh / year
            lcoeunit(j,2)=lcoeunit(j-1,2)+costs(idx(j),9)*rmb2us; % cost of connection to national grid
            lcoeunit(j,3)=lcoeunit(j-1,3)+costs(idx(j),6)*rmb2us;% cost of substation
            lcoeunit(j,4)=lcoeunit(j-1,4)+costs(idx(j),7)*rmb2us; % cost of expanding power unit
            lcoeunit(j,5)=lcoeunit(j-1,5)+capacity*moduleprice2;   % cost of PV module
            lcoeunit(j,7)=lcoeunit(j-1,7)+costs(idx(j),12)*rmb2us;   % cost of land
            lcoeunit(j,8)=lcoeunit(j-1,8)-electricity*fossilfuel_emissionfactor*CO2_C; % abated annual CO2 emission from fossil fuel Mton C / year
            lcoeunit(j,9)=lcoeunit(j-1,9)-costs(idx(j),5)*1e-6;  % annual land carbon sink Mton C / year
            lcoeunit(j,10)=lcoeunit(j-1,10)+costs(idx(j),4)*1e-6;   % total land use change emission Mton C
        else
            CPPP(j,1)=capacity; % capacity potential MW
            SRR(j,1)=wind(idx(j),1);
            SRR(j,2)=wind(idx(j),2);            
            lcoeunit(j,1)=electricity; % electricity TWh / year
            lcoeunit(j,2)=costs(idx(j),9)*rmb2us; % cost of connection to national grid
            lcoeunit(j,3)=costs(idx(j),6)*rmb2us; % cost of substation
            lcoeunit(j,4)=costs(idx(j),7)*rmb2us; % cost of expanding power unit
            lcoeunit(j,5)=capacity*moduleprice2;   % cost of PV module
            lcoeunit(j,7)=costs(idx(j),12)*rmb2us;   % cost of land
            lcoeunit(j,8)=-electricity*fossilfuel_emissionfactor*CO2_C; % abated annual CO2 emission from fossil fuel Mton C / year
            lcoeunit(j,9)=-costs(idx(j),5)*1e-6;  % emission from reduced annual land carbon sink Mton C / year
            lcoeunit(j,10)=costs(idx(j),4)*1e-6;   % total land use change emission Mton C
        end
        lcoeunit(j,11)=lcoeunit(j,2)*(1+OMratio_majorline*discount1yr); % cost of connection to national grid
        lcoeunit(j,12)=lcoeunit(j,3)*(1+OMratio_substation*discount1yr); % cost of substation
        lcoeunit(j,13)=lcoeunit(j,4)*(1+OMratio_substation*discount1yr); % cost of distance in power unit
        lcoeunit(j,14)=lcoeunit(j,5)*(1+OMratio_substation*discount1yr); % cost of PV module
        lcoeunit(j,16)=lcoeunit(j,7)*(1+OMratio_substation*discount1yr); % cost of land
        % LCoE
        lcoeunit(j,20)=sum(lcoeunit(j,11:19),2)/lcoeunit(j,1)/degrat1yr/1000; % LCoE million USD2019/TWh->USD2019/kWh        
        lcoeunit1(j,1)=sum(lcoeunit(j,12:19),2)/(lcoeunit(j,1))/degrat1yr/1000; % LCoE million USD2019/TWh->USD2019/kWh   
        if lcoeunit(j,20)<lcoemin
            jopt=j;
            lcoemin=lcoeunit(j,20);
        end
        if lcoeunit1(j,1)<lcoemin1 & lcoeunit1(j,1)>0
            jopt1=j;
            lcoemin1=lcoeunit1(j,1);
        end
    end
    % plot LCoE cost
    capacity_potential=sum(costs(idx,11),1); % MW
    j=jopt;
    powers(i,1:20)=lcoeunit(j,1:20); % 20 for LCoE USD2019/kWh
    powers(i,21:29)=lcoeunit(j,11:19)/lcoeunit(j,1)/degrat1yr/1000; % breakdown of LCoE USD2019/kWh
    powers(i,30)=sum(costs(idx(1):idx(j),11),1); % power capacity MW     
    CPPP111(i,1)=CPPP(j,1); % 20 for LCoE USD2019/kWh
    j1=jopt1;
    powers1(i,1)=lcoeunit1(j1,1); % 20 for LCoE USD2019/kWh
    windd(i,1)=SRR(j,1)/SRR(j,2);
end
lcoe_onshorewind=powers(:,20);
aass(:,1)=powers(:,30);
aass(:,2)=powers(:,20);
save('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\lcoe_onshorewind.mat','lcoe_onshorewind'); % USD2019/kWh

%%
rate_energy=0.02; % per year
load('G:\Code1123\powerdemand_pro2060_2.mat')  
E_120_2060 = powerdemand_pro2060_2*1000/365; % TWh/year-> GWh/grid/day
E_120 = E_120_2060./((1+rate_energy)^(2060-2019)); % GWh/grid/day 2019年
load('G:\Code1123\ID_Pro_CN.mat');  % ID_Pro_CN(4800x1950)
load('G:\Code1123\ID_County_CN.mat'); % county_CN(4800x1950)
land_sink=load('G:\Code1123\land_sink_sr2.txt'); % land_sink(72x46)
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

[B,IX]=sort(powers(:,20),1); 
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
    optpowerunit(i,1:30)=powers(i2,1:30); % lat lon
    optpowerunit(i,31)=i2; % id in powerunit
    optpowerunit(i,32)=optpowerunit(i,8)+optpowerunit(i,9);  % net CO2 emissions Mt C / year = (Abatement-Landsink)
    if i>1
        optpowerunit(i,33)=optpowerunit(i-1,33)+optpowerunit(i,32); % net CO2 emissions Mt C / year
        optpowerunit(i,34)=optpowerunit(i-1,34)+optpowerunit(i,1); % electricity TWh / year
    else
        optpowerunit(i,33)=emission2060+optpowerunit(i,32); % net CO2 emissions Mt C / year = FF+(Abatement-Landsink)
        optpowerunit(i,34)=optpowerunit(i,1); % electricity TWh / year
    end
end

optpowerunit_onshorewind = optpowerunit(:,1:34);
save('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\optpowerunit_onshorewind.mat','optpowerunit_onshorewind'); % power ID
powerunit_IX_onshorewind = powerunit_IX;
save('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit_IX_onshorewind.mat','powerunit_IX_onshorewind'); % power ID
save('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\tranmission_lines_IX.mat','lines_IX'); %
powerunit_num_IX_onshorewind = powerunit_num_IX;
save('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit_num_IX_onshorewind.mat','powerunit_num_IX_onshorewind'); %