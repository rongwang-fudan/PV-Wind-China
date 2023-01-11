tic;
clear
% 1.农；2.工业；3.建筑业；4.服务；5.其他; 6.交通；
% 7.sapce heating；8.cooling；9.cooking; 10.water heating; 11.rural household electric appliance
load('H:\China C neutrality\Data\powerdemand_ele_c.mat'); 
load('H:\China C neutrality\ANS\powergenerat_monhour_pv_8TW_pro.mat')  % 288*34 % TWh/h  
load('H:\China C neutrality\ANS\powergenerat_monhour_onshorewind_8TW_pro.mat')  % 288*34 % TWh/h  
load('H:\China C neutrality\ANS\powergenerat_monhour_offshorewind_8TW_pro.mat')  % 288*34 % TWh/h  
powergenerat_monhour = powergenerat_monhour_pv_8TW_pro+powergenerat_monhour_onshorewind_8TW_pro+powergenerat_monhour_offshorewind_8TW_pro;
clear powergenerat_monhour_pv_8TW_pro
clear powergenerat_monhour_onshorewind_8TW_pro
clear powergenerat_monhour_offshorewind_8TW_pro
ratio = powergenerat_monhour./sum(powergenerat_monhour);
ratio_c = sum(powergenerat_monhour,2)/sum(sum(powergenerat_monhour));
[m,n]=find(isnan(ratio));
aa=unique(n);
for i = 1:size(aa,1)
    ratio(:,aa(i))=ratio_c;
end

% 和各省发电量分布有关
% 1.农；2.工业；3.建筑业；4.服务；5.其他; 9.cooking; 10.water heating; 11.rural household electric appliance
powerdemand1 = sum(powerdemand_ele_c(:,[1:5,9:11]),2); % TWh/year
for i = 1:1:34
    powerdemand_monhour1(:,i) = powerdemand1(i,1)/8760*288*ratio(:,i); % TWh/h
    powerdemand_monhour_1(:,i) = powerdemand_ele_c(i,1)/8760*288*ratio(:,i); % TWh/h
    powerdemand_monhour_2(:,i) = powerdemand_ele_c(i,2)/8760*288*ratio(:,i); % TWh/h
    powerdemand_monhour_3(:,i) = powerdemand_ele_c(i,3)/8760*288*ratio(:,i); % TWh/h
    powerdemand_monhour_4(:,i) = powerdemand_ele_c(i,4)/8760*288*ratio(:,i); % TWh/h
    powerdemand_monhour_5(:,i) = powerdemand_ele_c(i,5)/8760*288*ratio(:,i); % TWh/h
    powerdemand_monhour_9(:,i) = powerdemand_ele_c(i,9)/8760*288*ratio(:,i); % TWh/h
    powerdemand_monhour_10(:,i) = powerdemand_ele_c(i,10)/8760*288*ratio(:,i); % TWh/h
    powerdemand_monhour_11(:,i) = powerdemand_ele_c(i,11)/8760*288*ratio(:,i); % TWh/h
end

% heating和cooling,和温度数据有关
% 7.sapce heating；8.cooling；
load('H:\China C neutrality\Data\cool_pro_r2.mat'); % 各省每小时cooling耗能占比
load('H:\China C neutrality\Data\heat_pro_r2.mat'); % 各省每小时heating耗能占比
heat_pro_r2(:,8) = heat_pro_r2(:,5); % 海南和广东的供暖变化曲线一样
cool_pro_r2(:,12) = cool_pro_r2(:,15); % 黑龙江和吉林的制冷变化曲线一样

mmon = [31 28 31 30 31 30 31 31 30 31 30 31]; % 2017年
ti=0;
for i = 1:1:12
    monthh(ti+1:ti+mmon(i)*24,1)=i;
    ti=ti+mmon(i)*24;
end
di = 0;
for i = 1:12
    md = mmon(i)*24;
    for j = 1:md
        day1(j,1) = ceil(j/24);
    end
    dayy(di+1:di+md,1) = day1;
    di = di+md; 
    clear md
    clear day1
end
for i = 1:365
    hourr(i*24-23:i*24,1) = [1:1:24];
end
a=0;
for i = 1:1:12
    a=a+mmon(i);
    mmons(i,1)=a;
end
hourrs(2:13,1)= mmons.*24;
hourrs(1,1)= 0;
cool_pro_r2(find(isnan(cool_pro_r2)==1)) = 0;
heat_pro_r2(find(isnan(heat_pro_r2)==1)) = 0;

cool_pro_r2_season = zeros(288,34);
heat_pro_r2_season = zeros(288,34);
    for month = 1:1:12
        [m,n] = find(monthh==month);
        dayy1=dayy(m,1);
        hourr1=hourr(m,1);
        for hour = 1:1:24
            [mm,nn] = find(hourr1==hour);
            mm = mm + hourrs(month);
            for ii = 1:1:size(mm,1)
                cool_pro_r2_season(hour+(month-1)*24,:) = cool_pro_r2(mm(ii),:)+cool_pro_r2_season(hour+(month-1)*24,:);
                heat_pro_r2_season(hour+(month-1)*24,:) = heat_pro_r2(mm(ii),:)+heat_pro_r2_season(hour+(month-1)*24,:);
            end
        end
    end
    
for i = 1:1:34
    powerdemand_monhour_cool(:,i) = powerdemand_ele_c(i,8)/8760*288.*cool_pro_r2_season(:,i);
    powerdemand_monhour_heat(:,i) = powerdemand_ele_c(i,7)/8760*288.*heat_pro_r2_season(:,i);
end
powerdemand_monhour_cool(find(isnan(powerdemand_monhour_cool)==1)) = 0;
powerdemand_monhour_heat(find(isnan(powerdemand_monhour_heat)==1)) = 0;

% trans,和车流量数据有关系
% 6.交通；
load('H:\China C neutrality\Data\transele_monhour288.mat'); % 
for i = 1:1:34
    powerdemand_monhour_trans(:,i) = powerdemand_ele_c(i,6)/8760*288.*transele_monhour288(:,i);
end

%合并
powerdemand_monhour2060ele_xz = powerdemand_monhour1 + powerdemand_monhour_cool + powerdemand_monhour_heat + powerdemand_monhour_trans;
save('H:\China C neutrality\ANS\powerdemand_monhour2060ele_xz.mat','powerdemand_monhour2060ele_xz')  % 288*34 % TWh/h  
