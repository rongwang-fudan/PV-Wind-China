tic
clear;
load('G:\China C neutrality\Data\SR_12to18_day.mat'); % SR_12to18_day(4800x1950)
SR_12to18_day(find(isnan(SR_12to18_day)==1))=0;
% SOCIAL-ECONOMY 1
load('G:\China C neutrality\Data\GDP_CN2015_120.mat'); % GDP_CN2015_120(4800x1950) ten thousand yuan
load('G:\China C neutrality\Data\POP_CN2020.mat'); % POP_CN2020(4800x1950) capita
load('G:\China C neutrality\Data\ID_County_CN.mat'); % county_CN(4800x1950)
load('G:\China C neutrality\Data\ID_Pro_CN.mat'); % ID_Pro_CN(4800x1950)
load('G:\China C neutrality\Data\costland_PV120.mat');  % yuan/grid

% ENERGY Production 2
load('G:\China C neutrality\Data\CP_PV_120_all.mat'); % CP_PV_120_all(4800x1950)
% Capacity power of PV KW
load('G:\China C neutrality\Data\Ph_12to18_day.mat'); % Ph_12to18_day(4800x1950)
% Electricity by PV as an average of 2012 and 2018 (GWh/grid/day)
% Electricity = Area * Ratio_Panel_Grid * Solar Radiation * Factor_shield *
% 16.19% * Factor_temperature * 80.56%
% load('..\data\wind_kineticenergy_100m_mean.mat'); % E_100m_mean(4800x1950)
% Wind kinetic energy for 100-m eight after accounting for slope 1979-2018 (kWh/yr)
% Density 4D*14D, D=126 meters
% 59.3% of kinetic energy is convertable to electricity
% t=(Ph_12to18_day./CP_PV_120_all); image(t,'cdatamapping','scaled'); % kWh/year/m2

% Energy Consumption 3
load('G:\China C neutrality\Data\electricity_demand_2019.mat'); % E_120(4800x1950) electricity demand in 2019
E_120=E_120.*(100/365); % 1e8 kWh / year -> GWh/grid/day
load('G:\China C neutrality\Data\Gasoline_120.mat'); % Gasoline_120(4800x1950)  kg/grid/year
Gasoline_120=Gasoline_120.*(1/0.7/5.5*20*1e-6/365); % (kg/grid/year) / (kg/L) / (L/100km) * (kWh/100km) -> GWh/grid/day
load('G:\China C neutrality\Data\Diesel_120.mat');  % Diesel_120(4800x1950)
Diesel_120=Diesel_120.*(1/0.8/5.5*20*1e-6/365); % (kg/grid/year) / (kg/L) / (L/100km) * (kWh/100km) -> GWh/grid/day
E_120=E_120+Gasoline_120+Diesel_120;

% HIGH-RESOLUTION 4
% load('..\bigdata\P_PV_monhour_all_mean.mat'); % LC_CN(4800x1950)

% GEOGRAPHIC DATA 5
load('G:\China C neutrality\Data\LC_CN.mat'); % LC_CN(4800x1950)
% land use (%) 1-5 forests ENF/EBF/DNF/DBF/MF; 6 Closed Shrubland
% 7 open Shrubland; 8 woody savanna; 9 savanna; 10 grassland; 11 wetland
% 12 cropland; 13 urban; 14 cropland / vegetation mosaics; 15 snow / ice
% 16 barren; 17 water bodies; 255 unknown
load('G:\China C neutrality\Data\Carbon_density12030.mat'); % Carbon_density12030(4800x1950)
% Lai L, et al. 2016 Carbon emissions. Science Advances t C/ha
land_sink=load('G:\China C neutrality\Data\land_sink_sr2.txt'); % land_sink(72x46)
% LC_CN(72x46) gC/m2/yr land carbon sink from atmospheric inversion
% Wang et al.Large Chinese land carbon sink estimated from atmospheric carbon dioxide data. Nature, 2021

% Transmission of electricity 5
load('G:\China C neutrality\PV_power potential\ANS_PV1\powerunit.dat','-mat');
load('G:\China C neutrality\PV_power potential\ANS_PV1\unitid.dat','-mat');
load('G:\China C neutrality\PV_power potential\ANS_PV1\unitid_all.dat','-mat');
load('G:\China C neutrality\PV_power potential\ANS_PV1\tranmission_lines.dat','-mat'); % lines
% image(unitid,'cdatamapping','scaled');

% Grid network 7
xnet=zeros(4800,1950);
ynet=zeros(4800,1950);
gridarea=zeros(4800,1950);
Rearth    =  6371.3;      % km average radium of the earth
E_120cy=zeros(2373,1); % electricity demand in a county GWh/day
for i=1:4800
    for j=1:1950
        xnet(i,j)=i;
        ynet(i,j)=j;
        gridarea(i,j)=abs(Rearth^2*(sin(((55-(i-1)/120+1/240)+1/120)*pi/180)-sin((55-(i-1)/120+1/240)*pi/180))*1/30*pi/180); %km2
        cy=county_CN(i,j);
        if cy>=1 && cy<=2373
            E_120cy(cy,1)=E_120cy(cy,1)+E_120(i,j); % GWh/day/county
        end
    end
end
% t=(E_120./gridarea).*100; image(t,'cdatamapping','scaled'); % kWh/year/m2
% t=(Ph_12to18_day./gridarea).*365;image(t,'cdatamapping','scaled'); % kWh/year/m2
% t=(costland_PV120./gridarea);image(t,'cdatamapping','scaled'); % yuan/km2

% Suitability Factor of Land Use 8
idx0=find(LC_CN==11 | LC_CN==17);
Ph_12to18_day(idx0)=0; % excluding water land
suitablity = 0.15*ones(4800,1950);
idx01=find(LC_CN==15);
suitablity(idx01)=0.1; % excluding water land

Ph_12to18_day=Ph_12to18_day.* suitablity;
CP_PV_120_all=CP_PV_120_all.* suitablity;

costs=zeros(100000,12);
SR2=zeros(100000000,2); % 1列是日照时长，2列是格子数
costunits=zeros(4800,1950); % id of the cost unit
id=0; % id of the cost unit
% For a power unit, we divide the grids into cost units based on the
% distance to the center of the power unit. In this way, we can find the
% cheapest grids for solar energy
linecost = 2.676; % million RMB cost of transmission line per 1 km
% linecost1000kV = 7.001; % million RMB cost of 1000kV transmission line per 1 km
trascapacity = 300; % capacity of transformer MW
trascost = trascapacity * 0.912; % million RMB for transformer
idxm=find(lines(:,6)==0); numlines=idxm(end); % id of lines
topon=[1 2 2 4 4];
iidd=0;
for i=1:size(powerunit,1)
    display(i);
    idx=find(unitid==i);
    idx_all=find(unitid_all==i);
    centerx=powerunit(i,1); % lat 1/120
    centery=powerunit(i,2); % lon 1/30
    cy=county_CN(centerx,centery);
    powercy=E_120cy(cy,1) / topon(powerunit(i,4)); % electricity demand divided by the number of power units in this county GWh/day/county
    isnottransmit=0;
    z=floor(abs(xnet(idx)-centerx)+abs(ynet(idx)-centery)*4)+1;
    z_all=floor(abs(xnet(idx_all)-centerx)+abs(ynet(idx_all)-centery)*4)+1;
    zmax=max(z,[],1); % the longest distance (1/120 degree)
    unitlength=sqrt(gridarea(centerx,centery)*4)/4; % length of a 1/120 deg km
    cumulpower=0; % cumulative power
    cumulpower2=0; % cumulative power
    transid=1;
    id1=id+1; % the first cost unit for this power unit
    % Initializing costs
    costs(id1,6)=trascost; % cost of substation million RMB
    for j=1:zmax
        iidd = iidd+1;
        id=id+1;
        costs(id,1)=i; % id of the power unit
        costs(id,2)=transid; % id of the transformer unit
        costs(id,7)=2*pi*linecost*j*unitlength; % cost of line million RMB
        costs(id,8)=j; % distance
        idx2=find(z==j);
        idx222=find(z_all==j);        
        if size(idx2,1)>0
            for k=1:size(idx222,1)
                x=xnet(idx_all(idx222(k)));
                y=ynet(idx_all(idx222(k)));
                costs(id,12)=costs(id,12)+costland_PV120(x,y)/1e6*suitablity(x,y); % cost of purchasing the land million RMB                     
            end
            for k=1:size(idx2,1)
                x=xnet(idx(idx2(k)));
                y=ynet(idx(idx2(k)));
                cumulpower=cumulpower+CP_PV_120_all(x,y)/1e3; % capacity (kW) -> MW        
                SR2(id,1) = SR2(id,1)+SR_12to18_day(x,y);
                SR2(id,2) = SR2(id,2)+1;
                if cumulpower>trascapacity
                    id=id+1;
                    transid=transid+1;
                    cumulpower=cumulpower-trascapacity;
                    costs(id,6)=trascost; % cost of substation million RMB
                    SR2(id,1) = SR2(id-1,1);
                    SR2(id,2) = SR2(id-1,2);                    
                end
                % connection to the major line
                if isnottransmit==0 % already connected to the major line
                    cumulpower2=cumulpower2+Ph_12to18_day(x,y); % GWh/grid/day
                    if cumulpower2>powercy
                        id=id+1;
                        costs(id,9)=lines(numlines+i,8)*unitlength*linecost; % cost of transmission to the major line million RMB
                        costs(id,10)=lines(numlines+i,5); % id of major line (1-47)
                        isnottransmit=1;
                        SR2(id,1) = SR2(id-1,1);
                        SR2(id,2) = SR2(id-1,2);                    
                    end
                end
                costs(id,1)=i; % id of the power unit
                costs(id,2)=transid; % id of the transformer unit
                costs(id,3)=costs(id,3)+Ph_12to18_day(x,y)*0.365; % electricity (GWh/grid/day) -> TWh / year  
                costs(id,8)=j; % distance                              
                costs(id,11)=costs(id,11)+CP_PV_120_all(x,y)/1e3; % capacity (kW) -> MW            
                % Carbon for vegetation only
                % For wind energy, there is no land-use change emission or land carbon sink change
                if LC_CN(x,y)<=10 || LC_CN(x,y)==12 || LC_CN(x,y)==14
                    costs(id,4)=costs(id,4)+Carbon_density12030(x,y)*gridarea(x,y)*100*suitablity(x,y); % land-use change emission t C/ha -> ton C
                    costs(id,5)=costs(id,5)+land_sink(floor((73+(y-1)/30+180)/5)+1,floor((55-(x-1)/120+90)/4)+1)*gridarea(x,y)*suitablity(x,y); % land carbon sink gC/m2/yr -> ton C/yr
                end
                costunits(x,y)=id; % id of cost units
            end
        end
    end
end
SR2=SR2(1:id,:); % delete unused data

t=costs;
costs=t(1:id,:); % delete unused data
save('G:\China C neutrality\PV_power potential\ANS_PV1\costunits01.dat','costs');
save('G:\China C neutrality\PV_power potential\ANS_PV1\SR_12to18_day_power.dat','SR2');

% clear
% toc



