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
%Parallel
%channel_flag(:,3)=1;
%channel_flag(:,5)=1;
%U
%channel_flag(1:x_m-6,3)=1;
%channel_flag(1:x_m-6,5)=1;
%channel_flag(x_m-6,3:5)=1;
%N
%channel_flag(1:x_m-6,2)=1;
%channel_flag(7:x_m,6)=1;
%channel_flag(x_m-6,3)=1;
%channel_flag(7,5)=1;
%Geo
channel_flag(1:7,2)=1;
channel_flag((x_m-6):x_m,2)=1;
channel_flag(7,3:5)=1;
channel_flag(x_m-6,3:5)=1;
channel_flag(7:(x_m-6),3)=1;
channel_flag(7:(x_m-6),5)=1;

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
velo(:,5)=-velo(:,5);
%velo(:,2)=2*velo(:,2);
%velo(7,:)=-2*velo(x_m-6,:);
%velo(x_m-6,:)=-2*velo(x_m-6,:);
mass_rate=abs(velo*density_water*a_cross);
Re = density_water*abs(velo_in)*d_h/viscocity_water;
if Re>2300
    Nu= 0.023*Re^0.8*Pr^0.4;
else
    Nu=6;
end
h_w = Nu*lamda_water/d_h;
%% Simulation
% iteration steps
it=10000;
for i = 1:it
    t_pv(1,1) = (I_m+h_a*t_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*(t_pv(1,2)*...
        x_m_d/y_m_d+t_pv(2,1)*y_m_d/x_m_d))/(h_a*(y_m_d*z_d+2*x_m_d*y_m_d)...
        +lamda_pv*z_d*(x_m_d/y_m_d+y_m_d/x_m_d));
    %t_pv(1,1) = (I_m+h_a*t_a*(y_m_d*z_d+x_m_d*y_m_d)+lamda_pv*z_d*(t_pv(1,2)*...
    %    x_m_d/y_m_d+t_pv(2,1)*y_m_d/x_m_d))/(h_a*(y_m_d*z_d+x_m_d*y_m_d)...
    %    +lamda_pv*z_d*(x_m_d/y_m_d+y_m_d/x_m_d));
    
    x=1;
    for y = 2:y_m-1
        t_pv(x,y) = (I_m+h_a*t_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*...
            ((t_pv(x,y+1)+t_pv(x,y-1))*x_m_d/y_m_d+t_pv(x+1,y)*y_m_d/x_m_d)...
            +channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
            /(h_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*...
            (2*x_m_d/y_m_d+y_m_d/x_m_d)+channel_flag(x,y)*h_w*x_m_d*p_cross);
        %t_pv(x,y) = (I_m+h_a*t_a*(y_m_d*z_d+x_m_d*y_m_d)+lamda_pv*z_d*...
        %    ((t_pv(x,y+1)+t_pv(x,y-1))*x_m_d/y_m_d+t_pv(x+1,y)*y_m_d/x_m_d)...
        %    +channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
        %    /(h_a*(y_m_d*z_d+x_m_d*y_m_d)+lamda_pv*z_d*...
        %    (2*x_m_d/y_m_d+y_m_d/x_m_d)+channel_flag(x,y)*h_w*x_m_d*p_cross);
        
        if channel_flag(x,y)==1
                if velo(x,y)<0
                    if channel_flag(x+1,y)==1
                        t_up = t_water(x+1,y);
                    elseif channel_flag(x,y+1)==1
                        t_up = t_water(x,y+1);
                    end
                    t_water(x,y)=(t_pv(x,y)*h_w*x_m_d*p_cross+cp_water*...
                        mass_rate(x,y)*t_up)/(cp_water*mass_rate(x,y)+...
                        h_w*x_m_d*p_cross);
                end
        end
    end
    
    t_pv(1, y_m) = (I_m+h_a*t_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*(t_pv(1,y_m-1)*...
        x_m_d/y_m_d+t_pv(2,y_m)*y_m_d/x_m_d))/(h_a*(y_m_d*z_d+2*x_m_d*y_m_d)...
        +lamda_pv*z_d*(x_m_d/y_m_d+y_m_d/x_m_d));
    %t_pv(1, y_m) = (I_m+h_a*t_a*(y_m_d*z_d+x_m_d*y_m_d)+lamda_pv*z_d*(t_pv(1,y_m-1)*...
    %    x_m_d/y_m_d+t_pv(2,y_m)*y_m_d/x_m_d))/(h_a*(y_m_d*z_d+x_m_d*y_m_d)...
    %    +lamda_pv*z_d*(x_m_d/y_m_d+y_m_d/x_m_d));
    
    for x = 2:x_m-1
        t_pv(x,1) = (I_m+h_a*t_a*(2*x_m_d*y_m_d)+lamda_pv*z_d*((t_pv(x-1,1)...
            +t_pv(x+1,1))*y_m_d/x_m_d+t_pv(x,2)*x_m_d/y_m_d))...
            /(h_a*(2*x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+2*y_m_d/x_m_d));
        t_pv(x,y_m) = (I_m+h_a*t_a*(2*x_m_d*y_m_d)+lamda_pv*z_d*((t_pv(x-1,y_m)...
            +t_pv(x+1,y_m))*y_m_d/x_m_d+t_pv(x,y_m-1)*x_m_d/y_m_d))...
            /(h_a*(2*x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+2*y_m_d/x_m_d));
        %t_pv(x,1) = (I_m+h_a*t_a*(x_m_d*y_m_d)+lamda_pv*z_d*((t_pv(x-1,1)...
        %    +t_pv(x+1,1))*y_m_d/x_m_d+t_pv(x,2)*x_m_d/y_m_d))...
        %    /(h_a*(x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+2*y_m_d/x_m_d));
        %t_pv(x,y_m) = (I_m+h_a*t_a*(x_m_d*y_m_d)+lamda_pv*z_d*((t_pv(x-1,y_m)...
        %    +t_pv(x+1,y_m))*y_m_d/x_m_d+t_pv(x,y_m-1)*x_m_d/y_m_d))...
        %    /(h_a*(x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+2*y_m_d/x_m_d));
    end
    
    t_pv(x_m,1) = (I_m+h_a*t_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*...
        (t_pv(x_m,2)*x_m_d/y_m_d+t_pv(x_m-1,1)*y_m_d/x_m_d))/...
        (h_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+y_m_d/x_m_d));
    t_pv(x_m,y_m) = (I_m+h_a*t_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*...
        (t_pv(x_m,y_m-1)*x_m_d/y_m_d+t_pv(x_m-1,y_m)*y_m_d/x_m_d))/...
        (h_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+y_m_d/x_m_d));
    
    x=x_m;
    for y = 2:y_m-1
        t_pv(x,y) = (I_m+h_a*t_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*...
            ((t_pv(x,y+1)+t_pv(x,y-1))*x_m_d/y_m_d+t_pv(x-1,y)...
            *y_m_d/x_m_d)+channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
            /(h_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*...
            (2*x_m_d/y_m_d+y_m_d/x_m_d)+channel_flag(x,y)*h_w*x_m_d*p_cross);
        if channel_flag(x,y)==1
                if velo(x,y)>0
                    if channel_flag(x-1,y)==1
                        t_up = t_water(x-1,y);
                    elseif channel_flag(x,y-1)==1
                        t_up = t_water(x,y-1);
                    end
                    t_water(x,y)=(t_pv(x,y)*h_w*x_m_d*p_cross+cp_water*...
                        mass_rate(x,y)*t_up)/(cp_water*mass_rate(x,y)+...
                        h_w*x_m_d*p_cross);
                end
        end
    end
    
    for x=2:x_m-1
        for y =2:y_m-1
            t_pv(x,y) = (I_m+h_a*t_a*2*x_m_d*y_m_d+lamda_pv*z_d*...
                ((t_pv(x,y-1)+t_pv(x,y+1))*x_m_d/y_m_d+...
                ((t_pv(x-1,y)+t_pv(x+1,y))*y_m_d/x_m_d))+...
                channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
                /(h_a*2*x_m_d*y_m_d+lamda_pv*z_d*2*(x_m_d/y_m_d+y_m_d/x_m_d)...
                +channel_flag(x,y)*h_w*x_m_d*p_cross);
            if channel_flag(x,y)==1
                if velo(x,y)>0
                    if channel_flag(x-1,y)==1
                        t_up = t_water(x-1,y);
                    elseif channel_flag(x,y-1)==1
                        t_up = t_water(x,y-1);
                    end
                elseif velo(x,y)<0
                    if channel_flag(x+1,y)==1
                        t_up = t_water(x+1,y);
                    elseif channel_flag(x,y+1)==1
                        t_up = t_water(x,y+1);
                    end
                end
                t_water(x,y)=(t_pv(x,y)*h_w*x_m_d*p_cross+cp_water*...
                    mass_rate(x,y)*t_up)/(cp_water*mass_rate(x,y)+...
                    h_w*x_m_d*p_cross);
            end
        
        end 
    end
    
    
end
%%Result
figure;
contourf(t_pv(6:(x_m-5),:));


