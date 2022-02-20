tic
clear;

trans = load('G:\China C neutrality\Data\transmissionline0626.txt');
lines=zeros(10000,15); % start_y; start_x; end_y; end_x; id of major line (1 to 47 for 2025)
id=0; % id of transmission line
for i=1:size(trans,1)
    idx=find(trans(i,:)~=0);
    num=(size(idx,2)-1)/2;
    if num>2
        % first part
        id=id+1;
        lines(id,1:2)=trans(i,2:3);
        lines(id,3:4)=trans(i,6:7);
        lines(id,5)=i;
        % centeral part
        for j=1:(num-3)
            id=id+1;
            lines(id,1:4)=trans(i,(j*2+4):(j*2+7));
            lines(id,5)=i;
        end       
        % last part        
        id=id+1;
        lines(id,1:2)=trans(i,idx(end-1):idx(end));
        lines(id,3:4)=trans(i,4:5);
        lines(id,5)=i; 
    else
        id=id+1;
        lines(id,1:4)=trans(i,2:5);
        lines(id,5)=i;
    end
end
% N55-N15, E73-E138 lat1/120; lon1/30
load('G:\China C neutrality\Data\powerunit_offshorewind.mat'); % lat-y, lon-x, 所属省份
lines(1:id,1)=floor((55-lines(1:id,1))*120)+1; % conversion of lat to y
lines(1:id,3)=floor((55-lines(1:id,3))*120)+1; % conversion of lat to y
lines(1:id,2)=floor((lines(1:id,2)-73)*30)+1; % conversion of lon to x
lines(1:id,4)=floor((lines(1:id,4)-73)*30)+1; % conversion of lon to x
powerunit_offshorewind(:,1)=floor((55-powerunit_offshorewind(:,1))*120)+1; % conversion of lat to y
powerunit_offshorewind(:,2)=floor((powerunit_offshorewind(:,2)-73)*30)+1; % conversion of lon to x

load('G:\China C neutrality\Data\ID_Pro_CN.mat');  % ID_Pro_CN(4800x1950)
load('G:\China C neutrality\Data\ID_County_CN.mat'); % county_CN(4800x1950)
regions = [2,3,2,5,6,6,7,6,3,1,6,4,1,1,4,2,2,4,0,3,5,5,7,5,2,2,3,3,6,5,7,7,2,7];
num=size(powerunit_offshorewind,1);
for i=1:num
    mindistance=10000;
    substation=0;
    for j=1:id
        distance=abs(powerunit_offshorewind(i,1)-lines(j,1))+abs(powerunit_offshorewind(i,2)-lines(j,2))*4;
        if distance<mindistance
            prov=ID_Pro_CN(lines(j,3),lines(j,4));
            if prov>=1 && prov<=34
                mindistance=distance;
                substation=j;
            end
        end
    end
    lines(id+i,1:2)=powerunit_offshorewind(i,1:2); % start point of the line
    lines(id+i,3:4)=lines(substation,1:2); % end point of the line
    lines(id+i,5)=lines(substation,5); % major line
    lines(id+i,6)=substation;
    lines(id+i,7)=i;
    lines(id+i,8)=mindistance;
end

t=lines;
lines=t(1:(id+i),:); % delete unused data
for j=1:id
    lines(j,9)=county_CN(lines(j,1),lines(j,2));
    lines(j,10)=ID_Pro_CN(lines(j,1),lines(j,2));
    if lines(j,10)>=1 && lines(j,10)<=34
        lines(j,11)=regions(lines(j,10));
    end
    lines(j,12)=county_CN(lines(j,3),lines(j,4));
    lines(j,13)=ID_Pro_CN(lines(j,3),lines(j,4));
    if lines(j,13)>=1 && lines(j,13)<=34
        lines(j,14)=regions(lines(j,13));
    end
end
pro_ID = [3 5 6 8 9 17 18 25 26 28 33];
for j=(id+1):(id+i)
    lines(j,9)=county_CN(lines(j,1),lines(j,2));
    lines(j,10)=pro_ID(j-id);
    if lines(j,10)>=1 && lines(j,10)<=34
        lines(j,11)=regions(lines(j,10));
    end
    lines(j,12)=county_CN(lines(j,3),lines(j,4));
    lines(j,13)=ID_Pro_CN(lines(j,3),lines(j,4));
    if lines(j,13)>=1 && lines(j,13)<=34
        lines(j,14)=regions(lines(j,13));
    end
end

Rearth    =  6371.3;      % km average radium of the earth
for j=1:id
    x=lines(j,3);
    gridarea=abs(Rearth^2*(sin(((55-(x-1)/120+1/240)+1/120)*pi/180)-sin((55-(x-1)/120+1/240)*pi/180))*1/30*pi/180); %km2
    unitlength=sqrt(gridarea*4)/4; % length of a 1/120 deg km
    lines(j,15)=unitlength * sqrt((lines(j,3)-lines(j,1))^2+16*(lines(j,4)-lines(j,2))^2); % line distance km
end
save('G:\China C neutrality\Offshore wind_power potential\ANS_OFFS1\tranmission_lines.dat','lines');

% clear
% toc


