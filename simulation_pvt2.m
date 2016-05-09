%% 2-d simulation for pvt
%Xingliang Fang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

%% Define Geometry and Mesh
% Overall dimensions
x_d=1.900;   %length in longtitute direction, m
y_d=0.071;   %length in latitute direction
z_d=0.002;   %length in thickness direction
% Meshing Demisions
y_m=7;      %Element number in x direction
y_m_d=y_d/y_m;  %Element length in x direction
x_m=round(x_d/y_m_d);        %Element number in y direction with similar element length
x_m_d=x_d/x_m; %Element length in y direction
% Define channel geometry
channel_flag=zeros(x_m,y_m);
channel_type='N';   %Channel type: 'P'-parallel, 'U', 'N'
draw_channel;
figure
subplot(1,2,1)
contourf(channel_flag,1)

d_h=2.106/1000; %hydraulic diameter, m
a_cross=10.882/1000000;  %cross section area, m^2
p_cross=20.6693/1000;
%% Material Properties
h_a = 10;
lamda_pv= 202.4;
lamda_water=0.677;
cp_water= 4200;
density_water= 1000;
viscocity_water=0.001002;
Pr = 7;
%% Initialization
t_pv_ini=25;    %initial temperature for pv
t_pv=t_pv_ini*ones(x_m,y_m);
t_water_ini=25; %initial temperature for water
t_water=t_water_ini*channel_flag;
%ambient temperature
t_a =25;        
%solar radiation
I=1000;
I_m=I*x_m_d*y_m_d;  %elemental I
%flow velocity
velo_in=0.1;  %inlet velocity in m/s
velo=velo_in*channel_flag;  %identical flow velocity distribution
switch channel_type
    case 'P'
        %
    case 'U'
        velo(:,[5,7])=-velo(:,[5,7]);
    case 'N'
        velo(:,4)=-velo(:,4);
end
mass_rate=abs(velo*density_water*a_cross);
Re = density_water*abs(velo_in)*d_h/viscocity_water;
if Re>2300
    Nu= 0.023*Re^0.8*Pr^0.4;
else
    Nu=6;
end
h_w = Nu*lamda_water/d_h;

%%simulation
% iteration steps
it=100;
simulation_body;

%%Result
subplot(1,2,2)
contourf(t_pv(6:(x_m-5),:))


