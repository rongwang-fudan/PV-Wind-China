%%
tic
clear;
xnet=zeros(4800,1950);
ynet=zeros(4800,1950);
for i=1:4800
    for j=1:1950
        xnet(i,j)=i;
        ynet(i,j)=j;
    end
end
load('G:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit.dat','-mat');
load('G:\China C neutrality\Onshore wind_power potential\ANS_ONS1\unitid.dat','-mat');
load('G:\China C neutrality\Onshore wind_power potential\ANS_ONS1\tranmission_lines.dat','-mat'); % lines
load('G:\China C neutrality\Onshore wind_power potential\ANS_ONS1\costunits_fig1b.dat','-mat');
numpowerunit=size(powerunit,1);
rmb2us=1/6.8967; % RMB to USD2019
moduleprice=2.7*rmb2us; % USD2019/Wp  Module
consprice=0.85*rmb2us; % USD2019/Wp 
gridprice=0.17*rmb2us; % USD2019/Wp
otherprice=0.11*rmb2us; % USD2019/Wp
moduleprice2= moduleprice+consprice+gridprice+ otherprice; %  USD2019/Wp 

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
OMratio_majorline=1/lifetime_power*1.2;
OMratio_substation=1/lifetime_power*1.2;
fossilfuel_emissionfactor=0.827; % tCO2/MWh Brander
CO2_C=0.2727;
carbonprice=0; % USD2019/tonCO2
powers=zeros(numpowerunit,30);
unitid_lcoe = zeros(4800,1950);

for i= 1:numpowerunit 
    idx_unitid=find(unitid==i);
    centerx=powerunit(i,1); % lat 1/120
    centery=powerunit(i,2); % lon 1/30
    z=floor(abs(xnet(idx_unitid)-centerx)+abs(ynet(idx_unitid)-centery)*4)+1;
    idx=find(costs(:,1)==i);
    lcoeunit=zeros(size(idx,1),20);
    lcoeunit_cz=zeros(size(idx,1),2);
    jopt=0;
    lcoemin=1000;
    for j=1:size(idx,1)
        electricity=costs(idx(j),3); % electricity TWh / year
        capacity=costs(idx(j),11); % capacity potential MW
        if j>1
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
        if j ==1
            lcoeunit_cz(j,2)=sum(lcoeunit(j,11:19),2)/lcoeunit(j,1)/degrat1yr/1000; % LCoE million USD2019/TWh->USD2019/kWh
            lcoeunit_cz(j,1)=lcoeunit(j,1);
        else
            lcoeunit_cz(j,2)=(sum(lcoeunit(j,11:19),2)-sum(lcoeunit(j-1,11:19),2))/(lcoeunit(j,1)-lcoeunit(j-1,1))/degrat1yr/1000; % LCoE million USD2019/TWh->USD2019/kWh
            lcoeunit_cz(j,1)=lcoeunit(j,1)-lcoeunit(j-1,1);
        end
        
        if lcoeunit(j,20)<lcoemin
            jopt=j;
            lcoemin=lcoeunit(j,20);
        end
    end
    capacity_potential=sum(costs(idx,11),1); % MW
    j=jopt;
    
    [m,n] = find(z<=j);
    idx_unitid_plant=idx_unitid(sub2ind(size(idx_unitid), m, n));
    unitid_lcoe(sub2ind(size(unitid_lcoe), xnet(idx_unitid_plant),ynet(idx_unitid_plant)))=unitid(sub2ind(size(unitid), xnet(idx_unitid_plant),ynet(idx_unitid_plant)));
        
    powers(i,1:20)=lcoeunit(j,1:20); % 20 for LCoE USD2019/kWh
    powers(i,21:29)=lcoeunit(j,11:19)/lcoeunit(j,1)/degrat1yr/1000; % breakdown of LCoE USD2019/kWh
    powers(i,30)=sum(costs(idx(1):idx(j),11),1); % power capacity MW     
    i
    clear idx_unitid_plant
end
%%
save('G:\China C neutrality\Onshore wind_power potential\ANS_ONS1\unitid_lcoe.dat','unitid_lcoe');
