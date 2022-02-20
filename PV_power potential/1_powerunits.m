tic
clear;
clear global;
global np0 npfraction

load('G:\China C neutrality\Data\ID_County_CN.mat'); % county_CN(4800x1950)
gridarea=zeros(4800,1950);
Rearth    =  6371.3;      % km average radium of the earth
% N55-N15, E73-E138 lat1/120; lon1/30
for i=1:4800
    for j=1:1950
        lat=55-(i-1)/120+1/240;
        lon=73+(j-1)/30-1/60; 
        gridarea(i,j)=abs(Rearth^2*(sin((lat+1/120)*pi/180)-sin(lat*pi/180))*1/30*pi/180); %km2
    end
end

% ENERGY 1
% load('..\data\I_12to18_day.mat'); % I_12to18_day(4800x1950)
% Solar Radiation as an average of 2012 and 2018  GWh/grid/day
load('G:\China C neutrality\Data\SR_12to18_day.mat'); % SR_12to18_day(4800x1950)
% Hours of solar radiation as an average of 2012 and 2018 (hour/day)
% load('..\data\Ph_12to18_day.mat'); % Ph_12to18_day(4800x1950)
% Electricity by PV as an average of 2012 and 2018 (GWh/grid/day)
% Electricity = Area * Ratio_Panel_Grid * Solar Radiation * Factor_shield *
% 16.19% * Factor_temperature * 80.56%
% load('..\data\wind_kineticenergy_100m_mean.mat'); % wind_kinetic_energy(4800x1950)
% Wind kinetic energy for 100-m eight after accounting for slope 1979-2018 (GWh/grid/day)
% Density 4D*14D, D=126 meters
% 59.3% of kinetic energy is convertable to electricity
% t=(I_12to18_day./gridarea).*365;image(t,'cdatamapping','scaled'); % kWh/year/m2
% t=(SR_12to18_day).*365;image(t,'cdatamapping','scaled'); % kWh/year/m2
% t=(wind_kinetic_energy./gridarea)/24*1000;image(t,'cdatamapping','scaled'); % W/m2

% GEOGRAPHIC DATA 2
load('G:\China C neutrality\Data\LC_CN.mat'); % LC_CN(4800x1950)
% land use (%) 1-5 forests ENF/EBF/DNF/DBF/MF; 6 Closed Shrubland
% 7 closed Shrubland; 8 woody savanna; 9 savanna; 10 grassland; 11 wetland
% 12 cropland; 13 urban; 14 cropland / vegetation mosaics; 15 snow / ice
% 16 barren; 17 water bodies; 255 unknown
load('G:\China C neutrality\Data\Slope_CN.mat'); % Slope_CN(4800x1950)
% Ground slope (%)
load('G:\China C neutrality\Data\T_2012-2018ave.mat'); % T_12to18(4800x1950)
% Temperature (C degree)
nature_reserve = imread('G:\China C neutrality\Data\nature_reserve_CN120.tif');
nature_reserve=double(nature_reserve);
species_reg = [21 22 23 24 25 26 35 41 46]; % Ecological reserve of species resources
for i = 1:size(species_reg,2)
    nature_reserve(nature_reserve==species_reg(i))=100;
end

% suitable land
suitableland=zeros(4800,1950);
% land use (%) 1-5 forests ENF/EBF/DNF/DBF/MF; 6 Closed Shrubland
% 7 Open Shrubland; 8 woody savanna; 9 savanna; 10 grassland; 11 wetland
% 12 cropland; 13 urban; 14 cropland / vegetation mosaics; 15 snow / ice
% 16 barren; 17 water bodies; 255 unknown
idx1=find(LC_CN==6 | LC_CN==7 | LC_CN==8 | LC_CN==9 | LC_CN==10 | LC_CN==15 | LC_CN==16);
s1=zeros(4800,1950); s1(idx1)=1;
idx2=find(s1==1 & Slope_CN<5 & T_12to18>0 & SR_12to18_day>4.2 & nature_reserve~=100);
suitableland(idx2)=1;
% image(suitableland,'cdatamapping','scaled');

% grid network
xnet=zeros(4800,1950);
ynet=zeros(4800,1950);
for i=1:4800
    for j=1:1950
        xnet(i,j)=i;
        ynet(i,j)=j;
    end
end

% number of grids as a function of distance (1/120 degree)
% np0=zeros(1000,1); np0(1)=1; 
% for i=1:999
%     np0(i+1)=np0(i)+i*4;
% end
idx=find(xnet<1800 & ynet<1800);
a=xnet(idx);
b=ynet(idx);
np0  =  npoint( a , b );

% Spatial Analysis
powerunit=zeros(10000,4);
id=0;
unitid=zeros(4800,1950); % suitable grid id of power unit in space
unitid_all=zeros(4800,1950); % all the grid id of power unit in space
npfraction=0.6;
for cy=1:2373
    % 5 types of topograph considered by R. Wang
    numgrid=zeros(5,13);
    numgrid_all=zeros(5,13);
    % TYPE 1
    idx=find(county_CN==cy & suitableland==1); %  
    idx2=find(county_CN==cy); %  
    if size(idx,1)<1
        nnn(cy) = 0;
        continue; % no grids
    end
    [ x1, y1, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
    numgrid(1,1)=num; % suitable grids number of the PV plant
    numgrid(1,2)=x1;
    numgrid(1,3)=y1;
    numgrid_all(1,1)=size(xy12,1);
    % Grid Information
    xy_suitabe=zeros(num,13*2);
    xy_suitabe(1:num,1:2)=xy; % suitable grids id of the PV plant
    xy_all=zeros(num,13*2);
    xy_all(1:size(xy12,1),1:2)=xy12; % All the grids id of the PV plant
    
    % TYPE 2
    for dom=1:2
        if dom==1
            idx=find(county_CN==cy & suitableland==1 & xnet>x1); %  Right
            idx2=find(county_CN==cy & xnet>x1); %  Right
        else
            idx=find(county_CN==cy & suitableland==1 & xnet<=x1); %  Left
            idx2=find(county_CN==cy & xnet<=x1); %  Left
        end
        if size(idx,1)<1
            continue; % no grids
        end
        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
        numgrid(2,1)=numgrid(2,1)+num;
        numgrid(2,9+dom)=num;
        numgrid(2,dom*2)=x;
        numgrid(2,dom*2+1)=y;
        numgrid_all(2,1)=numgrid_all(2,1)+size(xy12,1);
        numgrid_all(2,9+dom)=size(xy12,1);     
        xy_suitabe(1:num,(dom*2+1):(dom*2+2))=xy;
        xy_all(1:size(xy12,1),(dom*2+1):(dom*2+2))=xy12;
    end
    
    % TYPE 3
    for dom=1:2
        if dom==1
            idx=find(county_CN==cy & suitableland==1 & ynet>y1); %  Upper
            idx2=find(county_CN==cy & ynet>y1); %  Upper
        else
            idx=find(county_CN==cy & suitableland==1 & ynet<=y1); %  Bottom
            idx2=find(county_CN==cy & ynet<=y1); %  Bottom
        end
        if size(idx,1)<1
            continue; % no grids
        end
        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
        numgrid(3,1)=numgrid(3,1)+num;
        numgrid(3,9+dom)=num;
        numgrid(3,dom*2)=x;
        numgrid(3,dom*2+1)=y;
        numgrid_all(3,1)=numgrid_all(3,1)+size(xy12,1);
        numgrid_all(3,9+dom)=size(xy12,1);
        xy_suitabe(1:num,(dom*2+5):(dom*2+6))=xy;
        xy_all(1:size(xy12,1),(dom*2+5):(dom*2+6))=xy12;
    end
    
    % TYPE 4
    for dom=1:4
        if dom==1
            idx=find(county_CN==cy & suitableland==1 & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)>0); %  right-middle
            idx2=find(county_CN==cy & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)>0); %  right-middle
        elseif dom==2
            idx=find(county_CN==cy & suitableland==1 & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  down-middle
            idx2=find(county_CN==cy & (ynet-y1)<(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  down-middle
        elseif dom==3
            idx=find(county_CN==cy & suitableland==1 & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  left-middle
            idx2=find(county_CN==cy & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)<=0); %  left-middle
        elseif dom==4
            idx=find(county_CN==cy & suitableland==1 & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)>0); %  up-middle
            idx2=find(county_CN==cy & (ynet-y1)>=(xnet-x1) & (ynet-y1+xnet-x1)>0); %  up-middle
        end
        if size(idx,1)<1
            continue; % no grids
        end
        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
        numgrid(4,1)=numgrid(4,1)+num;
        numgrid(4,9+dom)=num;
        numgrid(4,dom*2)=x;
        numgrid(4,dom*2+1)=y;
        numgrid_all(4,1)=numgrid_all(4,1)+size(xy12,1);
        numgrid_all(4,9+dom)=size(xy12,1);
        xy_suitabe(1:num,(dom*2+9):(dom*2+10))=xy;
        xy_all(1:size(xy12,1),(dom*2+9):(dom*2+10))=xy12;
    end
    
    % TYPE 5
    for dom=1:4
        if dom==1
            idx=find(county_CN==cy & suitableland==1 & ynet>y1 & xnet<=x1); %  Upleft
            idx2=find(county_CN==cy & ynet>y1 & xnet<=x1); %  Upleft
        elseif dom==2
            idx=find(county_CN==cy & suitableland==1 & ynet>y1 & xnet>x1); %  Upright
            idx2=find(county_CN==cy & ynet>y1 & xnet>x1); %  Upright
        elseif dom==3
            idx=find(county_CN==cy & suitableland==1 & ynet<=y1 & xnet>x1); %  Downright
            idx2=find(county_CN==cy & ynet<=y1 & xnet>x1); %  Downright
        elseif dom==4
            idx=find(county_CN==cy & suitableland==1 & ynet<=y1 & xnet<=x1); %  Downleft
            idx2=find(county_CN==cy & ynet<=y1 & xnet<=x1); %  Downleft
        end
        if size(idx,1)<1
            continue; % no grids
        end
        [ x, y, num, xy, xy12 ]  =  fdist2( xnet(idx) , ynet(idx) , xnet(idx2) , ynet(idx2));  % Calculate the number of grids as a function of the distance to the center
        numgrid(5,1)=numgrid(5,1)+num;
        numgrid(5,9+dom)=num;
        numgrid(5,dom*2)=x;
        numgrid(5,dom*2+1)=y;
        numgrid_all(5,1)=numgrid_all(5,1)+size(xy12,1);
        numgrid_all(5,9+dom)=size(xy12,1);
        xy_suitabe(1:num,(dom*2+17):(dom*2+18))=xy;
        xy_all(1:size(xy12,1),(dom*2+17):(dom*2+18))=xy12;
    end
    
    % Find the best topography for power plants
    zmin=max(numgrid(:,1),[],1); idx1=find(numgrid(:,1)==zmin);
    if idx1(1)==1
        if numgrid(1,1)>2
            id=id+1; % 1 unit
            powerunit(id,1)=numgrid(1,2);
            powerunit(id,2)=numgrid(1,3);
            powerunit(id,3)=numgrid(1,1);
            powerunit(id,4)=1;
            for i=1:numgrid(1,1)
                unitid(xy_suitabe(i,1),xy_suitabe(i,2))=id;
            end
            for i=1:numgrid_all(1,1)
                unitid_all(xy_all(i,1),xy_all(i,2))=id;
            end
        end
    elseif idx1(1)==2
        for dom=1:2
            if numgrid(2,9+dom)>2
                id=id+1; % 2 unit
                powerunit(id,1)=numgrid(2,dom*2);
                powerunit(id,2)=numgrid(2,dom*2+1);
                powerunit(id,3)=numgrid(2,9+dom);
                powerunit(id,4)=2;
                for i=1:numgrid(2,9+dom)
                    unitid(xy_suitabe(i,dom*2+1),xy_suitabe(i,dom*2+2))=id;
                end
                for i=1:numgrid_all(2,9+dom)
                    unitid_all(xy_all(i,dom*2+1),xy_all(i,dom*2+2))=id;
                end
            end
        end
    elseif idx1(1)==3
        for dom=1:2
            if numgrid(3,9+dom)>2
                id=id+1; % 2 unit
                powerunit(id,1)=numgrid(3,dom*2);
                powerunit(id,2)=numgrid(3,dom*2+1);
                powerunit(id,3)=numgrid(3,9+dom);
                powerunit(id,4)=3;
                for i=1:numgrid(3,9+dom)
                    unitid(xy_suitabe(i,dom*2+5),xy_suitabe(i,dom*2+6))=id;
                end
                for i=1:numgrid_all(3,9+dom)
                    unitid_all(xy_all(i,dom*2+5),xy_all(i,dom*2+6))=id;
                end                
            end
        end
    elseif idx1(1)==4
        for dom=1:4
            if numgrid(4,9+dom)>2
                id=id+1; % 4 unit
                powerunit(id,1)=numgrid(4,dom*2);
                powerunit(id,2)=numgrid(4,dom*2+1);
                powerunit(id,3)=numgrid(4,9+dom);
                powerunit(id,4)=4;
                for i=1:numgrid(4,9+dom)
                    unitid(xy_suitabe(i,dom*2+9),xy_suitabe(i,dom*2+10))=id;
                end
                for i=1:numgrid_all(4,9+dom)
                    unitid_all(xy_all(i,dom*2+9),xy_all(i,dom*2+10))=id;
                end                
            end
        end
    elseif idx1(1)==5
        for dom=1:4
            if numgrid(5,9+dom)>2
                id=id+1; % 4 unit
                powerunit(id,1)=numgrid(5,dom*2);
                powerunit(id,2)=numgrid(5,dom*2+1);
                powerunit(id,3)=numgrid(5,9+dom);
                powerunit(id,4)=5;
                for i=1:numgrid(5,9+dom)
                    unitid(xy_suitabe(i,dom*2+17),xy_suitabe(i,dom*2+18))=id;
                end
                for i=1:numgrid_all(5,9+dom)
                    unitid_all(xy_all(i,dom*2+17),xy_all(i,dom*2+18))=id;
                end                
            end
        end
    end
    cy
end

t=powerunit;
powerunit=t(1:id,:); % delete unused data

save('G:\China C neutrality\PV_power potential\ANS_PV1\powerunit.dat','powerunit');
save('G:\China C neutrality\PV_power potential\ANS_PV1\unitid.dat','unitid');
save('G:\China C neutrality\PV_power potential\ANS_PV1\unitid_all.dat','unitid_all');

