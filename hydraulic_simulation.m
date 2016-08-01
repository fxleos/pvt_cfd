%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Numerical modelling and design of PVT modules               %%
%                      Semester Project                                   %
%                       Xingliang Fang                                    %
%             Chair of Architecture and Building Systems                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hydraulic simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
u=zeros(x_m,y_m);   %Velocity x-component
v=zeros(x_m,y_m);   %Velocity y-component
m=zeros(x_m,y_m);   %Velocity module

%% Simulation
%Define inlet and outlet conditions
switch channel_type
    case 'P'
    %Parallel
        start_p=[1,3;1,5];      %Inlet cell
        end_p=[x_m,3;x_m,5];    %Outline cell
        m_in=[1;1];             %Inlet velocity module
        velo_in=[1,0;1,0];      %Inlet velocity
    case 'P141'
    %Parallel 1-4-1
    start_p=[1,7;7,6;7,8;10,5;10,10;9,4;9,11;x_m-8,6;x_m-8,9;x_m-5,8];
    end_p=[7,7;9,5;9,10;x_m-8,5;x_m-8,10;x_m-8,5;x_m-8,10;x_m-6,8;x_m-6,8;x_m,8];
    m_in=[1;0.5;0.5;0.25;0.25;0.25;0.25;0.5;0.5;1];
    velo_in=[1,0;0,-1;0,1;1,0;1,0;0,-1;0,1;0,1;0,-1;1,0];
    case 'U'
    %U
    start_p=[1,1];
    end_p=[x_m,1];
    m_in=[1];
    velo_in=[1,0];
    case 'N'
    %N
    start_p=[1,1];
    end_p=[x_m,y_m];
    m_in=[1];
    velo_in=[1,0];
    case 'N12N'
    %N12
    start_p=[1,1];
    end_p=[1,y_m];
    m_in=[1];
    velo_in=[1,0];
    case 'MEBU'
    %MEBU
    start_p=[1,6];
    end_p=[2,2];
    m_in=[1];
    velo_in=[1,0];
end
%Calculating velocity field
[x_in,y_in]=size(start_p);
if x_in==1
    [u,v,m]=velocity_water( velo_in,m_in,channel_flag,start_p,end_p );
else
    for i=1:x_in
        [u_t,v_t,m_t]=velocity_water( velo_in(i,:),m_in(i),channel_flag,start_p(i,:),end_p(i,:) );
        u=u+u_t;
        v=v+v_t;
        m=m+m_t;
    end
end

%% Post-calculating mass flow rate and Reynolds number in each element
mass_rate=m*velo_amp*density_water*a_cross;
Re = density_water*(abs(u)+abs(v))*velo_amp*d_h/viscocity_water;