%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Numerical modelling and design of PVT modules               %%
%                      Semester Project                                   %
%                       Xingliang Fang                                    %
%             Chair of Architecture and Building Systems                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

%% Define Geometry and Mesh

%Overall dimensions
x_d=1.900;      %length in longtitute direction, m
y_d=0.071*2;    %length in latitute direction,m
z_d=0.002;      %length in thickness direction, m

%Select channel type from geometry pool
channel_type='P141';    %Channel type: 'P','U','N','P141','N12N','MEBU'
%Mesh the geometry
meshing;
figure
subplot(2,1,1)
contourf(channel_flag',1)

%Channel cross-section
d_h=2.106/1000;             %hydraulic diameter, m
a_cross=10.882/1000000;     %cross-sectional area, m^2
p_cross=20.6693/1000;       %cross-sectional perimeter,m

%% Thermal and hydarulic Properties
lamda_pv= 202.4;            %Conductivity of PV modules
lamda_water=0.677;          %Conductivity of working fluid
cp_water= 4200;             %Specific heat of working fluid
density_water= 1000;        %Density of working fluid
viscocity_water=0.001002;   %Viscocity of working fluid
Pr = 7;                     %Prandtl number of working fluid

%% Electric Parameters

%NOCT conditions
Area_noct=1.956*0.992;      %PV Module area
T_noct=25;                  %NOCT temperaure
I_noct=1000;                %Nominal Solar irradiation
V_op_noct=46;               %Nominal open circuit voltage
I_sc_noct=9.17;             %Nominal short circuit current
P_max_noct=325;             %Nominal maximum power output

%Temperature coefficients
c_Vop=-0.29/100;
c_Isc=0.05/100;
c_Pmax=-0.39/100;

%PV Type
pv_type='MONO';         %'MONO' or 'CIGS'

%PV placement
pv_placement='long';    %'long' or 'lati'

%PV connection
pv_connect='s';         %'s'-series, 'p'-parallel

%% Boundary conditions

%Equivalent heat transfer coefficient with ambient air
h_a = 10;
%Ambient air temperature
t_a =25;

%Solar irradiation
I=1000;
I_m=I*x_m_d*y_m_d;  %Calculated elemental I

%Inlet flow velocity
velo_amp=0.1;
%Inlet fluid temperature
t_water_ini=25;

%% Simulation parameters
%initial temperature for pv
t_pv_ini=25;
%Iteration steps
it=1000;
