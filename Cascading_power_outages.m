% Installation requirement: The code is running in MATLAB R2022a.
% Install Gurobi 10.1, YALMIP R20230609, MATPOWER 7.0.

% The methodology for simulating the wind field of Hurricane Fiona in 2022 draws upon the foundational work of Chavas, Lin, and Emanuel (2015), who developed a comprehensive model to depict the full radial structure of the tropical cyclone wind field. This model provides a physics-based method to simulate hurricane wind profiles. The reference for this method is as follows:
% Chavas, D., N. Lin, and K. Emanuel (2015). A model for the complete radial structure of the tropical cyclone wind field. Part I: Comparison with observed structure. J. Atmos. Sci. 

clear,clc

gurobi_setup();

if ~exist('yalmip', 'file')
    error('YALMIP is required');
end

% Region index
               % 1        2            3       4           5      6       7 
Region_Name = ["Ponce", "Mayaguez","Caguas","Arecibo","Bayamon","Carolina", "SanJuan"];



Simu_times = 1000;

% Roughness 
Rough_table = readtable('Windtable05mean.csv');
lat_R = Rough_table.LAT;
lon_R = Rough_table.LON;
Roughness_R = Rough_table.MEAN;
points = [lat_R,lon_R];

% Using KDTreeSearcher
kdtree = KDTreeSearcher(points);


% Hazard intensity

[sizex,sizey] = size(Wind{1});
size_t = length(Wind);

height0 = 10; % Average wind speed height
height = 10; % distribution network
height_trans = 30.48; % 100-foot towers
height_solar = 3; %


Wind_bus0 = zeros(size_t,length(LAT_bus));
Wind_bus = zeros(size_t,length(LAT_bus));
Rain_bus = zeros(size_t,length(LAT_bus));

Wind_bus0_trans = zeros(size_t,length(LAT_bus));
Wind_bus_trans = zeros(size_t,length(LAT_bus));

Wind_bus_solar = zeros(size_t,length(LAT_bus));

% %mpc = loadcase(PuertoRicoGrid)  % grid infrastructure

for k = 1:length(LAT_bus)


        lat = LAT_bus(k);
        lon = LON_bus(k);

%         [~,index_lat] = (min(abs(lat_R - lat)));
%         [~,index_lon] = (min(abs(lon_R - lon)));
% 
%         [row_r, col_r] = worldToIntrinsic(RR, x_R(index_lat), y_R(index_lon));
%         %roughness(k) = RasterRough(row_r, col_r);


        idx = knnsearch(kdtree, [lat, lon], 'K', 1);
        
        nearest_point_info = Rough_table(idx, :);

        roughness(k) = Roughness_R(idx);



    for t = 1:size_t


        % Wind calculate
        [Wind_bus0(t,k), Wind_bus(t,k), wind_direct, Radius] = wind_speed_cal(lat, lon, GLAT{t}, GLON{t}, Wind{t}, roughness(k), height0, height,Wind_U{t},Wind_V{t});
        [Wind_bus0_trans(t,k), Wind_bus_trans(t,k), wind_direct, Radius] = wind_speed_cal(lat, lon, GLAT{t}, GLON{t}, Wind{t}, roughness(k), height0, height_trans,Wind_U{t},Wind_V{t});
        [Wind_bus0_solar(t,k), Wind_bus_solar(t,k), wind_direct, Radius] = wind_speed_cal(lat, lon, GLAT{t}, GLON{t}, Wind{t}, roughness(k), height0, height_solar,Wind_U{t},Wind_V{t});

    end

    a = 1;

end


% interpolate wind in every 10 min
time_slot = 6;
Timeframe = 1:24;

num0 = 1:time_slot:time_slot*length(Timeframe);
num_int = 1:time_slot*length(Timeframe)+1-time_slot;

for k = 1:length(LAT_bus)
    
    Wind_int0(:,k) = interp1(num0,Wind_bus0(:,k),num_int,'linear');
    Wind_int(:,k) = interp1(num0,Wind_bus(:,k),num_int,'linear');
    Rain_int(:,k) = interp1(num0,Rain_bus(:,k),num_int,'linear');

    Wind_int0_trans(:,k) = interp1(num0,Wind_bus0_trans(:,k),num_int,'linear');
    Wind_int_trans(:,k) = interp1(num0,Wind_bus_trans(:,k),num_int,'linear');

    Wind_int0_solar(:,k) = interp1(num0,Wind_bus0_solar(:,k),num_int,'linear');
    Wind_int_solar(:,k) = interp1(num0,Wind_bus_solar(:,k),num_int,'linear');
    
    
end


%% Radiation Decay Factor

LAT_center = Fiona_t(:,5);
LON_center = Fiona_t(:,6);

LAT_center_int = interp1(num0,LAT_center,num_int,'linear');
LON_center_int = interp1(num0,LON_center,num_int,'linear');

%parameter from Ceferino et al. (2023)
R0_int = interp1(num0,cell2mat(R0),num_int,'linear')/1000;  % unit in km
ROCI = 0.18.*R0_int + 226;

Speed_ib = Fiona_t(:,7).* 1.15078; %knot to mph
Speed_ib = interp1(num0,Speed_ib,num_int,'linear');
for i = 1:length(Speed_ib)

    if Speed_ib(i)<=95
        Category(i) = 1;
    elseif Speed_ib(i)<=110
        Category(i) = 2;l
    elseif Speed_ib(i)<=129
        Category(i) = 3;
    elseif Speed_ib(i)<=156
        Category(i) = 4;
    else
        Category(i) = 5;
    end
    
end


%parameter from Ceferino et al. (2023)
for k = 1:length(LAT_bus)
    for t = 1:length(LAT_center_int)
                
    [arclen_grid,az_grid] = distance(LAT_bus(k),LON_bus(k),LAT_center_int(t),LON_center_int(t));

    distance_bus(t,k) = deg2km(arclen_grid);
    R_dist(t,k) = distance_bus(t,k)/ ROCI(t);

    if (R_dist(t,k) + (-0.126*Category(t)+1.15)) / (2.48-0.139*Category(t)) <=1
        factor_decay(t,k) = (0.0965*Category(t)+ 1.97) * log( (R_dist(t,k) + (-0.126*Category(t)+1.15)) / (2.48-0.139*Category(t))   );
    else
    end

    exp_decay(t,k) = exp(factor_decay(t,k));

    end
end

%% ms to mph

% m/s to mph
%----------------------------------------------%
Wind_bus_mph = Wind_bus.*2.23693629;
wind_bus_gustmph = Wind_bus_mph.*gust_f;
Wind_gust = Wind_int.*2.23693629.*gust_f;

Wind_gust_trans = Wind_int_trans.*2.23693629.*gust_f;
Wind_gust_solar = Wind_int_solar.*2.23693629.*gust_f;

% without roughness
Wind_bus_mph0 = Wind_bus0.*2.23693629;
wind_bus_gustmph0 = Wind_bus_mph0.*gust_f;
Wind_gust0 = Wind_int0.*2.23693629.*gust_f;

Wind_gust0_trans = Wind_int0_trans.*2.23693629.*gust_f;
Wind_gust0_solar = Wind_int0_solar.*2.23693629.*gust_f;


%----------------------------------------------%



% Calcualte Maxmium/Mean wind at each bus
%----------------------------------------------%

for t = 1:length(num_int)
    for i = 1: length(Wind_gust0(1,:) )

        Max_wind_bus(t,i) = max(Wind_gust(1:t,i)); 
        Max_wind_bus_trans(t,i) = max(Wind_gust_trans(1:t,i)); 

        Max_wind_bus0(t,i) = max(Wind_gust0_trans(1:t,i));

        Max_wind_bus_solar(t,i) = max(Wind_gust_solar(1:t,i));

        
    end
end




%% Vulnerability  

standard_230 = 300;
standard_115 = 200;

num_tower = zeros(length(dist_line),1);

for i = 1:length(voltage_line)

    if voltage_line(i) == 230
        
        num_tower(i) = dist_line(i) .*1000 ./ standard_230;
        
    else
        num_tower(i) = dist_line(i) .*1000 ./ standard_115;

    end

end

num_tower = ceil(num_tower);
num_line = num_tower - 1;

num_tower_from = ceil(num_tower/2);
num_tower_to = num_tower - num_tower_from;
num_line_from = ceil(num_line/2);
num_line_to = num_line - num_line_from;


%% fragility function

Fragility = readtable('Fragility_function.csv');
Fragility_wind_ms = Fragility.WindSpeed_ms_1_;
Fragility_wind_mph = Fragility_wind_ms(91:178).* 2.23693629 

Prob_line = Fragility.TransmissionLines(91:178);

fcn_line = @(b, x) cdf('Lognormal', x, b(1), b(2));
[B_line, resnorm_line] = fminsearch(@(b) norm(Prob_line - fcn_line(b, Fragility_wind_mph)), [1; 1]);
evf_line = fcn_line(B_line, Fragility_wind_mph);

x_line = linspace(min(Fragility_wind_mph), 300, 2000);
evf_fine = fcn_line(B_line, x_line);

% Plot the results
plot(Fragility_wind_mph, Prob_line, '.');
hold on;
plot(x_line, evf_fine, '-r');
hold off;
grid;
xlabel('3-s Wind Gust (m/s)');
ylabel('Failure percentage');
title('Combined Lognormal Fit for All Regions');

% Calculate R^2
y_mean = mean(Prob_line);
SStot = sum((Prob_line - y_mean).^2);
SSres = sum((Prob_line- evf_line).^2);
R2_line = 1 - SSres/SStot;
disp(['R^2 = ', num2str(R2)]);

text(median(Fragility_wind_mph), median(Prob_line)/2, ...
    sprintf('\\mu = %6.3f\n\\sigma = %6.3f\nR^2 = %6.3f', B_line(1), B_line(2), R2_line));



%Solar radiation decay - utility-sacle solar
                            %coto   horizon  AES     Oriana  San Fermin
index_utility_solar =         [56,     62,      64,     69     74];
index_utility_solar_connect = [ 4,      8,      9 ,     21,    37];

% start with 8pm UTC (20:00 - next day 20:00), local time UTC-4
solar_percent = [ 0 0 0 0 0 0 0 0 0 0 10 20 40 70 90 100 100 100 90 80 60 40 20 10];

gen_solar  = mpc.gen(index_utility_solar-51,:);
gen_solar_output  = mpc.gen(index_utility_solar-51,2);

cap_solar = gen_solar(:,9);


ratepower_utility = 0.5;
num_panel_soalr = ceil(cap_solar.*1000./ ratepower_utility);


solar_percent_int = interp1(num0,solar_percent,num_int,'linear');  % unit in km

for t = 1:length(solar_percent_int)
    gen_solar_output_int(t,:) = gen_solar_output.*solar_percent_int(t)/100;

    gen_solar_output_decay(t,:) =  gen_solar_output_int(t,:).* exp_decay(t,index_utility_solar_connect);
end



%% 
% Read Feeder Data
feeder_data = cell(1, numel(Region_Name));

for i = 1:numel(Region_Name)
    sheet_name = Region_Name(i);
    feeder_data{i} = readtable('PRFeeder_LUMA0529_withLATLON.xlsx', 'Sheet', sheet_name);

    Demand_feeder{i} = feeder_data{i}.FeederPeakDemand_MVA;
    Total_Demand_feeder{i} = sum(Demand_feeder{i})
    
end


% Gust wind for each distribution feeder
for i = 1:7

    bus_index{i} = feeder_data{i}.BusID ;
    Wind_feeder{i} = Wind_gust(:,bus_index{i});

    Num_feeder{i} = length(bus_index{i});

    Max_Wind_feeder{i} = zeros(size(Wind_feeder{i}));
    %Max_Wind_feeder(:,0) = 
    for t = 1: length(Wind_feeder{i}(:,1))

        for k = 1:Num_feeder{i}

            Max_Wind_feeder{i}(t,k) = max(Wind_feeder{i}(1:t,k));
        

        end
    end

end


%%



% Population
People = [212508	216605	250514	177369	216496	141663	253068];
Total_People = sum(People);

% Initialize
Percent_Live_feeder = cell(Simu_times, 7);
Total_Percent_Live_feeder = cell(Simu_times, 1);


deltat = 1/6;
t2 = 23;
t1 = 0;
n = (t2-t1)/deltat
%n = time_slot*(length(Timeframe)-1)  

results = struct();


for num = 1:Simu_times % 1000 realizations

    for i = 1:7

        %Resistance_feeder = 0; %inverse function of the fragility curve

        random_feeder{i} = rand(1,Num_feeder{i}); %feeder

        Resistance_feeder_log{i} = icdf('Lognormal', random_feeder{i}, B{i}(1), B{i}(2));

        %exp
        y_inverse = -log(1-random_feeder{i});
        Resistance_feeder_exp{i} = (log(y_inverse / exp(b1(i))) / a1(i));


        %Compare max wind and resistance
        % using bsxfun to compare resilience/hazard intensity and store in binary_matrix 
        Failure_State_log{i} = double(bsxfun(@gt, Max_Wind_feeder{i}, Resistance_feeder_log{i}));
        Failure_State_exp{i} = double(bsxfun(@gt, Max_Wind_feeder{i}, Resistance_feeder_exp{i}));

        for k = 1: Num_feeder{i}
            
            Loss_Distribution_log{i}(:,k) = Failure_State_log{i}(:,k).*Demand_feeder{i}(k);
            Loss_Distribution_exp{i}(:,k) = Failure_State_exp{i}(:,k).*Demand_feeder{i}(k);

        end

        Total_Loss_feeder{i} = sum(Loss_Distribution_log{i}, 2);
        Percent_Loss_feeder{i} = Total_Loss_feeder{i}/Total_Demand_feeder{i};
        Percent_Live_feeder{num,i} = 100- Percent_Loss_feeder{i}*100;

        Live_customer(:,i) = Percent_Live_feeder{num,i}*People(i);

        %Output
        Out_Region_Percent{i}(:,num) = 100- Percent_Loss_feeder{i}*100;



    end

    Total_Live_customer = sum(Live_customer,2);  
    Total_Percent_Live_feeder{num} = Total_Live_customer./Total_People;

    %Output 
    Out_Total_Percent(:,num) = Total_Percent_Live_feeder{num}; 





    %% Resistance-based damage estimation
    Demand_bus = zeros(51,1);
    for i = 1:7

        uniqiue_bus_index{i} = unique(bus_index{i});
        %selected_bus_index = cell()
        for j = 1:length(uniqiue_bus_index{i} )

            col_bus_feeder = find(bus_index{i} == uniqiue_bus_index{i}(j));

            Demand_bus(uniqiue_bus_index{i}(j)) = sum(Demand_feeder{i}(col_bus_feeder));


            Loss_bus_feeder =  Loss_Distribution_log{i}(:,col_bus_feeder);
            Demand_loss_bus(:,uniqiue_bus_index{i}(j)) = sum(Loss_bus_feeder,2);
            Percent_demand_loss_bus(:,uniqiue_bus_index{i}(j)) = Demand_loss_bus(:,uniqiue_bus_index{i}(j))./Demand_bus(uniqiue_bus_index{i}(j));

            
        end

    end


    Percent_demand_live_bus = ones(size( Percent_demand_loss_bus)) - Percent_demand_loss_bus;


    %  Resistance values for tower and tranmsision line

    mpc00 = loadcase(PuertoRicoGrid)
    mpc0 = mpc00;
    line_from = mpc0.branch(:,1);
    line_to = mpc0.branch(:,2);

    syms w_t
    func_tower = 1 - exp (- alpha* exp(delta * w_t ) );



    Rmin_tower = zeros(1,length(dist_line));

    for i = 1:length(dist_line)
        %%-------------------tower

        y_tower_from = unifrnd(0,1,length(num_tower_from),1);
        y_tower_to = unifrnd(0,1,length(num_tower_to),1);

        Resist_tower_from = inv_func_tower(y_tower_from);
        Resist_tower_to = inv_func_tower(y_tower_to);

        Rmin_tower_from(i) = min(Resist_tower_from);
        Rmin_tower_to(i) = min(Resist_tower_to);

        %%-------------------line

        y_line_from = unifrnd(0,1,length(num_line_from),1);
        y_line_to = unifrnd(0,1,length(num_line_to),1);

        Resist_line_from = icdf('Lognormal', y_line_from, B_line(1), B_line(2));
        Resist_line_to = icdf('Lognormal', y_line_to, B_line(1), B_line(2));
        %---------------------%
        Rmin_line_from(i) = min(Resist_line_from);
        Rmin_line_to(i) = min(Resist_line_to);

    end

    index_from_bus = line_from(1:length(dist_line));
    index_to_bus = line_to(1:length(dist_line));

    Max_wind_bus_from = Max_wind_bus0(:,index_from_bus);
    Max_wind_bus_to = Max_wind_bus0(:,index_to_bus);
    
    %Compare max wind and resistance
    %%-------------------tower
    Failure_State_tower_from = double(bsxfun(@gt, Max_wind_bus_from,  Rmin_tower_from ));
    Failure_State_tower_to   = double(bsxfun(@gt, Max_wind_bus_to,  Rmin_tower_to ));

    %Failure_State_tower = Failure_State_tower_from + Failure_State_tower_to;
    Ione = ones(size(Failure_State_tower_from));

    Failure_State_tower =  Ione - ( Ione - Failure_State_tower_from ).* ( Ione - Failure_State_tower_to );
    
    
    
    %%-------------------line

    
    Failure_State_line_from = double(bsxfun(@gt, Max_wind_bus_from,  Rmin_line_from ));
    Failure_State_line_to   = double(bsxfun(@gt, Max_wind_bus_to,  Rmin_line_to ));

    
    Failure_State_line =  Ione - ( Ione - Failure_State_line_from ).* ( Ione - Failure_State_line_to );
    

    Failure_State_corridor =  Ione - ( Ione - Failure_State_line ).* ( Ione - Failure_State_tower );

    
 
    underground_line = find(underground_state ==1);
    % underground line will not fail
    for jj = underground_line'
        Failure_State_corridor(:,jj) = zeros(length(Failure_State_corridor),1);
    end

    state_corridor = Ione - Failure_State_corridor;



    %% Solar fragility - Utility
    % solar parameter from Ceferino et al. (2023)
    miu_utility = log(90);
    sigma_utility = 0.15;

    % Utility_Scale Solar Fragility+Decay
    for j = 1:length(index_utility_solar)

        Max_wind_solar{j} = Max_wind_bus0(:,index_utility_solar_connect(j));

        Failure_rate_utility_solar(:,j) = logncdf(Max_wind_solar{j}./2.23693629,miu_utility, sigma_utility);
        I_solar = ones(length(Failure_rate_utility_solar(:,j)),1);
        gen_solar_output_final(:,j) =  gen_solar_output_decay(:,j).*(I_solar-Failure_rate_utility_solar(:,j));

    end
    



%%   Solar radiation decay - distriubted rooftop solar  

    solar_percent_total = sum(solar_percent);
    
    solar_percent_average = mean(solar_percent);
    
    adjust_factor = capacity_factor/solar_percent_average
    
    solar_percent_adjust = adjust_factor.*solar_percent *0.01 % percent to decimal
    solar_percent_adjust_int = interp1(num0,solar_percent_adjust,num_int,'linear')
    
    for i = 1:numel(Region_Name)
        sheet_name = Region_Name(i);
        feeder_data{i} = readtable('PRFeeder_LUMA0529_withLATLON.xlsx', 'Sheet', sheet_name);
        % This part involves Critical Energy/Electric Infrastructure Information (CEII)
        % Request should be sent to LUMA Energy or Federal Energy Regulatory Commission 
    
        unique_sub = unique(feeder_data{i}.SubstationID,'rows');
        unique_bus = unique(feeder_data{i}.BusID);
    
        for k = 1:length(unique_bus)
    
            index_DG = (feeder_data{i}.BusID == unique_bus(k));
            sub_DG = sum(index_DG.*feeder_data{i}.ExistingDGCapacity_kW)/1000
            sub_Demand = sum(index_DG.*feeder_data{i}.FeederPeakDemand_MVA)
    
            Penetration(unique_bus(k),:) = [unique_bus(k), (sub_DG/sub_Demand)*100 , sub_DG, sub_Demand ]
    
        end
    
    end
    DG_capacity_bus =   Penetration(:,3)'   %139*1  * 1*51
    P_DG_normal = solar_percent_adjust_int' * DG_capacity_bus 
    P_DG_decay  = P_DG_normal .* exp_decay
    
    
    miu_roof = log(80);
    sigma_roof = 0.32;
    
    for j = 1: length(DG_capacity_bus )
    
        Max_wind_solar_rooftop = Max_wind_bus0(:,j) 
        Max_wind_solar_rooftop = Max_wind_bus_solar(:,j) 
        Failure_rate_roof_solar(:,j) = logncdf(Max_wind_solar_rooftop./2.23693629,miu_roof, sigma_roof);
        I_solar = ones(length(Failure_rate_roof_solar(:,j)),1);
       
    end
    
    Ieye = ones(size(Failure_rate_roof_solar))
    P_DG_decay_failure  = P_DG_decay .* (Ieye - Failure_rate_roof_solar)
    
    Load_increase_solar = P_DG_normal - P_DG_decay_failure;



    %% Cascading power outages
    % Initialization
    mpc(1) = mpc0;

    for t = 1:n+1
        mpc(t) = mpc0;
        % Check transmission line state
        mpc(t).branch(1:length(dist_line),11) = state_corridor(t,:)';  % 11th row is status
    end



    load_total0 = sum(mpc0.bus(:,3));
    load('Load_Profile0523.mat');
    load_percent = interp1(num0,Load_Profile(:,3),num_int,'linear')';  % unit in km
    load_total_t = load_percent.*load_total0;
 
    for t = 1:n
 
        L_base0(t,:) = mpc(67).bus(:,3)'.*load_percent(t);

        L_base(t,:) = L_base0(t,:);
        L_base(t,1:length(Load_increase_solar(t,:))) = L_base0(t,1:length(Load_increase_solar(t,:)) )  + Load_increase_solar(t,:);

    end


    bus_length = length(Percent_demand_live_bus(1,:))

    mpc(1).gen(:,2) = mpc(1).gen(:,2).*load_percent(1);
    mpc(1).bus(1:bus_length,3) = mpc0.bus(1:bus_length,3).*Percent_demand_live_bus(1,:)'.*load_percent(1);


    total_gen(1) = sum( mpc(1).gen(:,2))
    total_load(1) = sum(mpc(1).bus(:,3))

    % utility-scale solar geneartion
    gen_solar_index =  find( ismember(mpc(1).gen(:,1), index_utility_solar))
    for t = 1:n+1
        mpc(t).gen(gen_solar_index,2) = gen_solar_output_final(t,:)'
    end


    % Spatiotemporal power outages

    for t = 1:n  % n=138

         if t ==1 
                 L_c(t,:) = mpc(1).bus(1:bus_length,3)';
                 Cut(t,:) = zeros(size(L_c(t,:)));
         else

             for k = 1:bus_length

                 Cut(t,k) = 0;

                 if Percent_demand_live_bus(t,k) ==0 || L_base(t,k) ==0  || L_c(t-1,k)<= 0.01 

                     L_c(t,k) = 0;
                 else
                     L_c(t,k) = L_c(t-1,k)* (1 - Cut(t-1,k)/L_c(t-1,k) ) * ( Percent_demand_live_bus(t,k)/ Percent_demand_live_bus(t-1,k) ) ...
                                           * (L_base(t,:) / L_base(t-1,:) );

                 end


                 
             end
                
         end

         mpc(t).bus(1:bus_length,3) = L_c(t,:)';

        %----- Network Reconfiguration ---------%
        [mpc1_final,flag_failure(t),Condition_total(t,:),Condition_sum(t,:)] = NetworkReconfig(mpc(t));
        %-------------------------------------- %

        [groups_new,isolated_new] = find_islands(mpc1_final);




        if isempty(groups_new)
            total_load(t:n) = 0;
            total_gen(t:n) = 0 ;
            
            loss_percent(t:n) = 0;

            break
            
        else

            mpc1_sub = extract_islands(mpc1_final);


            total_load_sub = 0;
            total_gen_sub = 0;

            N_islands = size(groups_new,2)

            for i = 1:N_islands
                total_load_sub(i) = sum(mpc1_sub{1,i}.bus(:,3));
                total_gen_sub(i) = sum(mpc1_sub{1,i}.gen(:,2));
            end
            total_load(t,num) = sum(total_load_sub)
            total_gen(t,num) = sum(total_gen_sub)

        end

        Cut_load =  mpc(t).bus(1:bus_length,3)  - mpc1_final.bus(1:bus_length,3);
        Cut(t,:) = Cut_load';


        mpc(t+1).bus = mpc1_final.bus;


        mpc(t+1).gen = mpc1_final.gen;
        mpc(t+1).branch(:,11) = mpc(t+1).branch(:,11).*mpc1_final.branch(:,11); %corridor status




        loss_percent(t) = sum(L_c(t,:))/ sum(L_base(t,:))

     end


    A_result{num} = {};

    simulation = struct();
    simulation.loss_percent = loss_percent
    simulation.mpc = mpc
    simulation.Condition = Condition_sum
    simulation.failure_tower = Failure_State_tower
    simulation.failure_line = Failure_State_line


    A_result{num} = simulation

end

