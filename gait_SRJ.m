%% gait.m -- Starter File for Gait Data Analysis
%% B17 Biomechanics -- Hilary Term 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

% Load data into arrays Xmk (markers) and Xfp (forceplate) for Subject X
% You must replace 'X' below with the code for your Subject (A, B, C, D)
mk = csvread('Subject-D-markers.csv');
fp = csvread('Subject-D-forceplate.csv');

% Find the number of rows and columns in the input files
[rmk,cmk] = size(mk);    [rfp,cfp] = size(fp);

% initialize arrays of datapt number in marker file & assign to an array
% note that ':' means 'all the rows (or columns)'
datapt = zeros(rmk,1);
datapt(:,1) = mk(:,1);

% Initialize arrays for marker data (right leg only)
sac = zeros(rmk,3);      % sacrum (SACR)
R_asi = zeros(rmk,3);    % anterior superior iliac spine (LASI)
thi = zeros(rmk,3);      % thigh wand marker (THI) 
kne = zeros(rmk,3);      % lateral femoral condyle (KNE)
tib = zeros(rmk,3);      % tibia wand marker (TIB)
ank = zeros(rmk,3);      % later malleolus (LMA)
hee = zeros(rmk,3);      % heel (HEE)
toe = zeros(rmk,3);      % 2nd metatarsal head (TOE)

% Assign xyz coordinates of markers to right side sacrum, asis, thigh, knee,
% ankle, heel, and toe arrays
sac(:,1:3) = mk(:,2:4);
R_asi(:,1:3) = mk(:,8:10);             
thi(:,1:3) = mk(:,29:31);      
kne(:,1:3) = mk(:,32:34); 
tib(:,1:3) = mk(:,35:37);
ank(:,1:3) = mk(:,38:40);      
hee(:,1:3) = mk(:,41:43);      
toe(:,1:3) = mk(:,44:46);    

%%%%%%%%%%%%% YOU NEED TO CONTINUE THE CODE FROM HERE 

%% Question 1: Plot yz trajectories
     figure(1)
     hold on
     
     plot(sac(:,2),sac(:,3))
     text(sac(rmk,2),sac(rmk,3),'SACRUM')
     
     plot(R_asi(:,2),R_asi(:,3))
     text(R_asi(rmk,2),R_asi(rmk,3),'ILIAC')
     
     plot(thi(:,2),thi(:,3))
     text(thi(rmk,2),thi(rmk,3),'THIGH')
     
     plot(kne(:,2),kne(:,3))
     text(kne(rmk,2),kne(rmk,3),'KNEE')
     
     plot(tib(:,2),tib(:,3))
     text(tib(rmk,2),tib(rmk,3),'TIBIA')
     
     plot(ank(:,2),ank(:,3))
     text(ank(rmk,2),ank(rmk,3),'ANKLE')
     
     plot(hee(:,2),hee(:,3))
     text(hee(rmk,2),hee(rmk,3),'HEEL')
     
     plot(toe(:,2),toe(:,3))
     text(toe(rmk,2),toe(rmk,3),'TOE')
     
     legend('SACRUM','ILIAC','THIGH','KNEE','TIBIA','ANKLE','HEEL','TOE')
     xlabel('Y (Posterior - Anterior)')
     ylabel('Z (Inferior - Superior)')
     axis('equal')
     title('Subject D Marker Trajectories')


%% Question 2: Calculate instantaneous fowrard velocity of body

L_asi = zeros(rmk,3);        % initialise array for LASI marker
L_asi(:,1:3) = mk(:,5:7);    % assign xyz coordinates of marker

com = (R_asi + L_asi)/2;     % calculate array for centre of mass
com_y = com(:,2);            % extract y coordinates for centre of mass 

sample_rate = 100;           % define sample rate
dt = 1/sample_rate;          % define time interval

velocity = zeros(rmk-1,1);     % initialise velocity vector 

% Eulers method to calculate velocity in m/s
for k = 1:rmk-1

    velocity(k) = velocity(k) + (com_y(k+1) - com_y(k))/(dt*1000);  

end

% Calculate average velocity 
av_velocity = zeros(rmk-1,1);              % initialise average velocity vector
av_velocity(1:end) = mean(velocity);       % calculate average velocity 


% Plot forward velocity against data point number 
figure(2)
hold on
plot(datapt(1:end-1), velocity)
plot(datapt(1:end-1), av_velocity)
xlabel('Data Point')
ylabel('Forward Velocity (m/s)')
title('Forward Velocity of Centre of Mass')
legend('Forward Velocity','Average Velocity')



%% Question 3: Shank length

% Calculating shank length in 2D using Euclidean distance
shank_length_2D = zeros(rmk,1);        % initialise shank length 2D vector 

for k = 1:rmk 
    
    shank_length_2D(k) = shank_length_2D(k) + pdist([kne(k,2:3); ank(k,2:3)],'euclidean');
    
end

% Calculating shank length in 3D using Euclidean distance
shank_length_3D = zeros(rmk,1);     % initialise shank length 3D vector

for k = 1:rmk
    
    shank_length_3D(k) = shank_length_3D(k) + pdist([kne(k,:); ank(k,:)],'euclidean');
    
end
    
figure(3)
hold on
plot(datapt, shank_length_2D)
plot(datapt, shank_length_3D)
xlabel('Data Point')
ylabel('Shank Length (mm)')
title('Shank Length Variation')
legend('2D length','3D length')


%% Question 4: Define the gait cycle

% Assign variables within forceplate array
Fx = fp(:,2);
Fy = fp(:,3);
Fz = fp(:,4);
Mx = fp(:,5);
My = fp(:,6);
Mz = fp(:,7);

datapt2 = zeros(rfp, 1);    % initialise array to contain data points
datapt2(:,1) = fp(:,1);     % fill array with data points from force plate csv

% plot vertical components of marker trajectories of ankle, heel and toe
figure(4)
hold on
plot(datapt, ank(:,3))
plot(datapt, hee(:,3))
plot(datapt, toe(:,3))
xlabel('Data Point')
ylabel('Height (mm)')
title('Vertical Components of Marker Trajectories')
legend('Ankle','Heel','Toe')

% Plot vertical component of ground reaction force 
figure(5)
plot(datapt2, Fz)
xlabel('Data Point')
ylabel('Force')
title('Vertical Component of Ground Reaction Force')


% Finding specific data points in force plate data
index_nonzero = find(Fz);              % array of index numbers for all non-zero Forces
index_HS = min(index_nonzero);         % index number for first heel strike   
index_TO = 1 + max(index_nonzero);     % index number for first toe-off

DP_HS_1 = datapt2(index_HS);           % data point number for first heel strike
DP_TO = datapt2(index_TO);             % data point number for first toe-off

% Data points in marker data
DP_HS_2 = 460;                         % data point number for second heel strike (recorded from graph)


% Time for one gait cycle 
sample_rate2 = 1000;              % sample rate for force plate data
dt2 = 1/1000;                     % time interval for force plate data

t_HS1 = DP_HS_1 * dt2;            % time of first heel strike
t_HS2 = DP_HS_2 * dt;             % time of second heel strike 

gait_time = t_HS2 - t_HS1;        % time period for one gait cycle


% Time for stance phase 
t_TO = DP_TO * dt2;               % time of first toe-off

stance_time = t_TO - t_HS1;       % stance phase time


% Time for swing phase 
swing_time = gait_time - stance_time;

% Stance phase (% gait cycle)
stance_percentage = (stance_time / gait_time)*100;

% Swing phase (% gait cycle)
swing_percentage = 100 - stance_percentage;




%% Question 5: Stride length & Candence

% Calculating stride length 
DP_HS_1_marker = fix(DP_HS_1 / 10);                % equivalent marker data point for first heel strike
DP_HS_2_marker = DP_HS_2;                          % marker data point for second heel strike

index_HS_1 = find(mk(:,1) == DP_HS_1_marker);      % index number for first heel strike  
index_HS_2 = find(mk(:,1) == DP_HS_2_marker);      % index number for second heel strike

stride_length = (hee(index_HS_2,2) - hee(index_HS_1,2))/1000;   % stride length in metres
    

%Calculating cadence
cadence = (mean(velocity) / stride_length) * 60;    % cadence in steps/min 



%% Question 6: Knee joint angle

% Using prep work notation to define lengths
a = sqrt((kne(:,2) - thi(:,2)).^2 + (kne(:,3) - thi(:,3)).^2);      
b = sqrt((tib(:,2) - kne(:,2)).^2 + (tib(:,3) - kne(:,3)).^2);
c = sqrt((tib(:,2) - thi(:,2)).^2 + (tib(:,3) - thi(:,3)).^2);

theta = acos((a.^2 + b.^2 - c.^2) ./ (2 .* a .* b));     % calculating knee joing angle using cosine rule

% Index number for toe-off
DP_TO_marker = fix(DP_TO / 10);                 % equivalent marker data point for toe-off
index_TO = find(mk(:,1) == DP_TO_marker);       % index number for first heel strike

% Plot knee joint angle against data number 
figure(6)
plot(datapt, theta)
set(gca, 'YDir','reverse')
xlabel('Data Point')
ylabel('Angle (radians)')
title('Plot of Knee Joint Angle')
text(datapt(index_HS_1), theta(index_HS_1), '\leftarrow HS')
text(datapt(index_HS_2), theta(index_HS_2), '\leftarrow HS')
text(datapt(index_TO), theta(index_TO), '\leftarrow TO')


