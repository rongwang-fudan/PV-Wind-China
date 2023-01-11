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
    if sum(costs(idx(1):idx(jopt),11),1)<=100000 % power capacity MW
        j=jopt;
    else if sum(costs(idx(1):idx(jopt),11),1)>100000
            for iiiiii = 1:jopt
                cp_1230(iiiiii,1) = sum(costs(idx(1):idx(iiiiii),11),1);
            end
            [~,Index123] = min(abs(cp_1230-100000));
            j=Index123;
            clear cp_1230
        end
    end
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
        optpowerunit(i,34)=optpowerunit(i-1,34)+optpowerunit(i,1); % electricity TWh / year
    else
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