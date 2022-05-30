tic
clear;
clear global;
global np0 npfraction

load('H:\China C neutrality\Data\ID_County_CN.mat'); % county_CN(4800x1950)
gridarea=zeros(4800,1950);
Rearth=6371.3;% km average radium of the earth
% N55-N15, E73-E138 lat1/120; lon1/30
for i=1:4800
    for j=1:1950
        lat=55-(i-1)/120+1/240;
        lon=73+(j-1)/30-1/60; 
        gridarea(i,j)=abs(Rearth^2*(sin((lat+1/120)*pi/180)-sin(lat*pi/180))*1/30*pi/180); %km2
    end
end

% ENERGY 1
load('H:\China C neutrality\Data\windpower_100m_12to18_day_2.mat'); % wind_kineticenergy(4800x1950) (GWh/grid/day) for 100-m height
% P_h=S×ρ×CF×UTI_coef×ARR_coef÷1000
% S: Available land area（m2）
% ρ=P_Wp/(7D*5D): power density for the current turbine model, W/m2；P_Wp=2500kW；D=103m: Rotor diameter；
% CF: the capacity factor for the curr ent wind turbine model and hub height in the cell；
% UTI_coef= 0.95: the utilization efficiency, due to technical failures and so on (literature values ranged from 0.94 to 0.98)；
% ARR_coef= 0.90: the array efficiency factor, due to wake effects in the wind turbine arrays. (literature values range from 0.7 to 0.925).

% GEOGRAPHIC DATA 2
load('H:\China C neutrality\Data\LC_CN.mat'); % LC_CN(4800x1950)
% land use (%) 1-5 forests ENF/EBF/DNF/DBF/MF; 6 Closed Shrubland
% 7 closed Shrubland; 8 woody savanna; 9 savanna; 10 grassland; 11 wetland
% 12 cropland; 13 urban; 14 cropland / vegetation mosaics; 15 snow / ice
% 16 barren; 17 water bodies; 255 unknown
load('H:\China C neutrality\Data\Slope_CN.mat'); % Slope_CN(4800x1950)
% Ground slope (%)
load('H:\China C neutrality\Data\CF_mean.mat') %单位 km2

% suitable land
suitableland=zeros(4800,1950);
% land use (%) 1-5 forests ENF/EBF/DNF/DBF/MF; 6 Closed Shrubland
% 7 closed Shrubland; 8 woody savanna; 9 savanna; 10 grassland; 11 wetland
% 12 cropland; 13 urban; 14 cropland / vegetation mosaics; 15 snow / ice
% 16 barren; 17 water bodies; 255 unknown

nature_reserve = imread('H:\China C neutrality\Data\nature_reserve_CN120.tif');
nature_reserve=double(nature_reserve);
species_reg = [21 22 23 24 25 26 35 41 46]; % Ecological reserve of species resources
for i = 1:size(species_reg,2)
    nature_reserve(nature_reserve==species_reg(i))=100;
end
DEM = imread('H:\China C neutrality\Data\DEM_CN12030.tif');

idx1=find(LC_CN==6 | LC_CN==7 | LC_CN==8 | LC_CN==9 | LC_CN==10 | LC_CN==12 | LC_CN==14 | LC_CN==16);
s1=zeros(4800,1950); s1(idx1)=1;
idx2=find(s1==1 & Slope_CN<20 & CF>0.2 & nature_reserve~=100 & DEM<3000);
suitableland(idx2)=1;

% grid network
xnet=zeros(4800,1950);
ynet=zeros(4800,1950);
for i=1:4800
    for j=1:1950
        xnet(i,j)=i;
        ynet(i,j)=j;
    end
end

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

save('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\powerunit.dat','powerunit');
save('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\unitid.dat','unitid');
save('H:\China C neutrality\Onshore wind_power potential\ANS_ONS1\unitid_all.dat','unitid_all');

