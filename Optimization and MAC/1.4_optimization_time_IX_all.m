tic
clear;
rmb2us=1/6.8967; % RMB to USD2019
fossilfuel_emissionfactor=0.783;
CO2_C=0.2727;
lifetime_power=25;
discount=0.05; % per year
discount1yr=0;
for t=1:lifetime_power
    discount1yr=discount1yr+1/(1+discount)^(t-1);
end
% PV:
OMratio_majorline_PV=0.01; %;
OMratio_substation_PV=0.01;
% wind
OMratio_majorline_wind=0.03;
OMratio_substation_wind=0.03;

load('H:\China C neutrality\PV_power potential\ANS_PV1\tranmission_lines.dat','-mat');  % lines
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
lines(numlines+1:end,:)=[];

load('H:\China C neutrality\PV_power potential\ANS_PV1\optpowerunit_PV.mat'); %
load('H:\China C neutrality\PV_power potential\ANS_PV1\powerunit_IX_PV.mat'); %
load('H:\China C neutrality\PV_power potential\ANS_PV1\tranmission_lines_IX.mat');  % lines_IX
lines_IX(size(optpowerunit_PV,1)+1:end,:)=[];
lines_IX_PV = lines_IX;
load('H:\China C neutrality\PV_power potential\ANS_PV1\powerunit_num_IX_PV.mat');  %
optpowerunit_PV(:,35) = 1;
optpowerunit_PV(:,40) = powerunit_IX_PV; % power plant type ID
load('H:\China C neutrality\ANS\unitmin_module_pv.mat')
optpowerunit_PV(:,41) = unitmin_pv; % Buiding year

load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\optpowerunit_onshorewind.mat'); %
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit_IX_onshorewind.mat'); %
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\tranmission_lines_IX.mat');  % lines_IX
lines_IX(size(optpowerunit_onshorewind,1)+1:end,:)=[];
lines_IX_onshorewind = lines_IX;
load('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit_num_IX_onshorewind.mat');  %
optpowerunit_onshorewind(:,35) = 2;
optpowerunit_onshorewind(:,40) = powerunit_IX_onshorewind; % power plant type ID
load('H:\China C neutrality\ANS\unitmin_module_ons.mat')
optpowerunit_onshorewind(:,41) = unitmin_ons; % Buiding year

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
optpowerunit_offshorewind(:,40) = powerunit_IX_offshorewind; % power plant type ID
load('H:\China C neutrality\ANS\unitmin_module_off.mat')
optpowerunit_offshorewind(:,41) = unitmin_off; % Buiding year

optpowerunit = [optpowerunit_PV;optpowerunit_onshorewind;optpowerunit_offshorewind];
lines_IX = [lines_IX_PV(:,1:15);lines_IX_onshorewind;lines_IX_offshorewind];
powerunit_num_IX = [powerunit_num_IX_PV;powerunit_num_IX_onshorewind;powerunit_num_IX_offshorewind];
[B,IX]=sort(optpowerunit(:,20),1);
numpowerunit = size(optpowerunit,1);
for i=1:numpowerunit
    i2=IX(i);
    powerunit_IX(i,1)=i2;
    optpowerunit_IX(i,1:41)=optpowerunit(i2,1:41); %  optpowerunit_IX(i2,1); electricity used by the county TWh / year
    lines_IX_IX(i,1:15)=lines_IX(i2,1:15); % lat lon
    powerunit_num_IX_IX(i,1:4)=powerunit_num_IX(i2,1:4); % lat lon
end
unitmin = optpowerunit_IX(:,41);
save('H:\China C neutrality\ANS\unitmin_module.mat','unitmin');