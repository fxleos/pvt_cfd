%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Numerical modelling and design of PVT modules               %%
%                      Semester Project                                   %
%                       Xingliang Fang                                    %
%             Chair of Architecture and Building Systems                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Thermal simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
t_pv=t_pv_ini*ones(x_m,y_m);
t_water=t_water_ini*channel_flag;

%% Simulation
for i = 1:it
    x=1;
    
    y=1;
    if Re(x,y)>2300
        if t_pv(x,y)>t_water(x,y)
            Nu= 0.023*Re(x,y)^0.8*Pr^0.4;
        else
            Nu= 0.023*Re(x,y)^0.8*Pr^0.3;
        end
    else
        if t_pv(x,y)>t_water(x,y)
            Nu= 5;
        else
            Nu= 4.5;
        end
    end
    h_w = Nu*lamda_water/d_h;
    t_pv(x,y) = (I_m+h_a*t_a*(y_m_d*z_d+x_m_d*y_m_d)+...
        lamda_pv*z_d*(t_pv(x,y+1)*x_m_d/y_m_d+t_pv(x+1,1)*y_m_d/x_m_d)+...
        channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
        /(h_a*(y_m_d*z_d+x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+y_m_d/x_m_d)...
        +channel_flag(x,y)*h_w*x_m_d*p_cross);


        
    for y = 2:y_m-1
        t_pv(x,y) = (I_m+h_a*t_a*(y_m_d*z_d+x_m_d*y_m_d)+lamda_pv*z_d*...
            ((t_pv(x,y+1)+t_pv(x,y-1))*x_m_d/y_m_d+t_pv(x+1,y)*y_m_d/x_m_d)...
            +channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
            /(h_a*(y_m_d*z_d+x_m_d*y_m_d)+lamda_pv*z_d*...
            (2*x_m_d/y_m_d+y_m_d/x_m_d)+channel_flag(x,y)*h_w*x_m_d*p_cross);

        
    if channel_flag(x,y)==1
    switch 3*u(x,y)+v(x,y)
        case 3
            %t_up = t_water(x-1,y);
        case -3
            t_up = t_water(x+1,y);
        case 1
            t_up = t_water(x,y-1);
        case -1
            t_up = t_water(x,y+1);
        case 4
            t_up = (t_water(x-1,y)+t_water(x,y-1))/2;
        case 2
            t_up = (t_water(x-1,y)+t_water(x,y+1))/2;
        case -2
            t_up = (t_water(x+1,y)+t_water(x,y))/2;
        case -4
            t_up = (t_water(x+1,y)+t_water(x,y+1))/2;
    end
    %t_water(x,y)=(t_pv(x,y)*h_w*x_m_d*p_cross+cp_water*...
        %mass_rate(x,y)*t_up)/(cp_water*mass_rate(x,y)+...
        %h_w*x_m_d*p_cross);
    end
    end
    
    y=y_m;
    t_pv(x,y) = (I_m+h_a*t_a*(y_m_d*z_d+x_m_d*y_m_d)+...
        lamda_pv*z_d*(t_pv(x,y-1)*x_m_d/y_m_d+t_pv(x+1,y)*y_m_d/x_m_d)+...
        channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
        /(h_a*(y_m_d*z_d+x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+y_m_d/x_m_d)...
        +channel_flag(x,y)*h_w*x_m_d*p_cross);

    if channel_flag(x,y)==1
    switch 3*u(x,y)+v(x,y)
        case 3
            t_up = t_water(x-1,y);
        case -3
            t_up = t_water(x+1,y);
        case 1
            t_up = t_water(x,y-1);
        case -1
            t_up = t_water(x,y+1);
        case 4
            t_up = (t_water(x-1,y)+t_water(x,y-1))/2;
        case 2
            t_up = (t_water(x-1,y)+t_water(x,y+1))/2;
        case -2
            t_up = (t_water(x+1,y)+t_water(x,y))/2;
        case -4
            t_up = (t_water(x+1,y)+t_water(x,y+1))/2;
    end
    if Re(x,y)>2300
        if t_pv(x,y)>t_water(x,y)
            Nu= 0.023*Re(x,y)^0.8*Pr^0.4;
        else
            Nu= 0.023*Re(x,y)^0.8*Pr^0.3;
        end
    else
        if t_pv(x,y)>t_water(x,y)
            Nu= 5;
        else
            Nu= 4.5;
        end
    end
    h_w = Nu*lamda_water/d_h;
    t_water(x,y)=(t_pv(x,y)*h_w*x_m_d*p_cross+cp_water*...
        mass_rate(x,y)*t_up)/(cp_water*mass_rate(x,y)+...
        h_w*x_m_d*p_cross);
    end
    
    for x = 2:x_m-1
        y=1;
        t_pv(x,y) = (I_m+h_a*t_a*(2*x_m_d*y_m_d)+lamda_pv*z_d*((t_pv(x-1,1)...
            +t_pv(x+1,1))*y_m_d/x_m_d+t_pv(x,y+1)*x_m_d/y_m_d)+...
            channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
            /(h_a*(2*x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+2*y_m_d/x_m_d)...
            +channel_flag(x,y)*h_w*x_m_d*p_cross);
        if channel_flag(x,y)==1
        switch 3*u(x,y)+v(x,y)
            case 3
                t_up = t_water(x-1,y);
            case -1
                t_up = t_water(x+1,y);
            case 1
                t_up = t_water(x,y-1);
            case -1
                t_up = t_water(x,y+1);
            case 4
                t_up = (t_water(x-1,y)+t_water(x,y-1))/2;
            case 2
                t_up = (t_water(x-1,y)+t_water(x,y+1))/2;
            case -2
                t_up = (t_water(x+1,y)+t_water(x,y))/2;
            case -4
                t_up = (t_water(x+1,y)+t_water(x,y+1))/2;
        end
        if Re(x,y)>2300
            if t_pv(x,y)>t_water(x,y)
                Nu= 0.023*Re(x,y)^0.8*Pr^0.4;
            else
                Nu= 0.023*Re(x,y)^0.8*Pr^0.3;
            end
        else
            if t_pv(x,y)>t_water(x,y)
                Nu= 5;
            else
                Nu= 4.5;
            end
        end
        h_w = Nu*lamda_water/d_h;
        t_water(x,y)=(t_pv(x,y)*h_w*x_m_d*p_cross+cp_water*...
            mass_rate(x,y)*t_up)/(cp_water*mass_rate(x,y)+...
            h_w*x_m_d*p_cross);
        end
        
        for y =2:y_m-1
            t_pv(x,y) = (I_m+h_a*t_a*2*x_m_d*y_m_d+lamda_pv*z_d*...
                ((t_pv(x,y-1)+t_pv(x,y+1))*x_m_d/y_m_d+...
                ((t_pv(x-1,y)+t_pv(x+1,y))*y_m_d/x_m_d))+...
                channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
                /(h_a*2*x_m_d*y_m_d+lamda_pv*z_d*2*(x_m_d/y_m_d+y_m_d/x_m_d)...
                +channel_flag(x,y)*h_w*x_m_d*p_cross);
            if channel_flag(x,y)==1
                switch 3*u(x,y)+v(x,y)
                    case 3
                        t_up = t_water(x-1,y);
                    case -3
                        t_up = t_water(x+1,y);
                    case 1
                        t_up = t_water(x,y-1);
                    case -1
                        t_up = t_water(x,y+1);
                    case 4
                        t_up = (t_water(x-1,y)+t_water(x,y-1))/2;
                    case 2
                        t_up = (t_water(x-1,y)+t_water(x,y+1))/2;
                    case -2
                        t_up = (t_water(x+1,y)+t_water(x,y))/2;
                    case -4
                        t_up = (t_water(x+1,y)+t_water(x,y+1))/2;
                end
                if Re(x,y)>2300
                    if t_pv(x,y)>t_water(x,y)
                        Nu= 0.023*Re(x,y)^0.8*Pr^0.4;
                    else
                        Nu= 0.023*Re(x,y)^0.8*Pr^0.3;
                    end
                else
                    if t_pv(x,y)>t_water(x,y)
                        Nu= 5;
                    else
                        Nu= 4.5;
                    end
                end
                h_w = Nu*lamda_water/d_h;
                t_water(x,y)=(t_pv(x,y)*h_w*x_m_d*p_cross+cp_water*...
                    mass_rate(x,y)*t_up)/(cp_water*mass_rate(x,y)+...
                    h_w*x_m_d*p_cross);
            end
        
        end 
        
        y=y_m;
        t_pv(x,y) = (I_m+h_a*t_a*(2*x_m_d*y_m_d)+lamda_pv*z_d*((t_pv(x-1,1)...
            +t_pv(x+1,1))*y_m_d/x_m_d+t_pv(x,y-1)*x_m_d/y_m_d)+...
            channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
            /(h_a*(2*x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+2*y_m_d/x_m_d)...
            +channel_flag(x,y)*h_w*x_m_d*p_cross);

        if channel_flag(x,y)==1
            switch 3*u(x,y)+v(x,y)
                case 3
                    t_up = t_water(x-1,y);
                case -3
                    t_up = t_water(x+1,y);
                case 1
                    t_up = t_water(x,y-1);
                case -1
                    t_up = t_water(x,y+1);
                case 4
                    t_up = (t_water(x-1,y)+t_water(x,y-1))/2;
                case 2
                    t_up = (t_water(x-1,y)+t_water(x,y+1))/2;
                case -2
                    t_up = (t_water(x+1,y)+t_water(x,y))/2;
                case -4
                    t_up = (t_water(x+1,y)+t_water(x,y+1))/2;
            end
            if Re(x,y)>2300
                if t_pv(x,y)>t_water(x,y)
                    Nu= 0.023*Re(x,y)^0.8*Pr^0.4;
                else
                    Nu= 0.023*Re(x,y)^0.8*Pr^0.3;
                end
            else
                if t_pv(x,y)>t_water(x,y)
                    Nu= 5;
                else
                    Nu= 4.5;
                end
            end
            h_w = Nu*lamda_water/d_h;
            t_water(x,y)=(t_pv(x,y)*h_w*x_m_d*p_cross+cp_water*...
                mass_rate(x,y)*t_up)/(cp_water*mass_rate(x,y)+...
                h_w*x_m_d*p_cross);
        end
    end
    
    x=x_m;
    
    y=1;
    t_pv(x,y) = (I_m+h_a*t_a*(y_m_d*z_d+x_m_d*y_m_d)+...
        lamda_pv*z_d*(t_pv(x-1,y)*y_m_d/x_m_d+t_pv(x,y+1)*x_m_d/y_m_d)+...
        channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
        /(h_a*(y_m_d*z_d+x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+y_m_d/x_m_d)...
        +channel_flag(x,y)*h_w*x_m_d*p_cross);
    if channel_flag(x,y)==1
    switch 3*u(x,y)+v(x,y)
        case 3
            t_up = t_water(x-1,y);
        case -3
            t_up = t_water(x+1,y);
        case 1
            t_up = t_water(x,y-1);
        case -1
            t_up = t_water(x,y+1);
        case 4
            t_up = (t_water(x-1,y)+t_water(x,y-1))/2;
        case 2
            t_up = (t_water(x-1,y)+t_water(x,y+1))/2;
        case -2
            t_up = (t_water(x+1,y)+t_water(x,y))/2;
        case -4
            t_up = (t_water(x+1,y)+t_water(x,y+1))/2;
    end
    if Re(x,y)>2300
        if t_pv(x,y)>t_water(x,y)
            Nu= 0.023*Re(x,y)^0.8*Pr^0.4;
        else
            Nu= 0.023*Re(x,y)^0.8*Pr^0.3;
        end
    else
        if t_pv(x,y)>t_water(x,y)
            Nu= 5;
        else
            Nu= 4.5;
        end
    end
    h_w = Nu*lamda_water/d_h;
    t_water(x,y)=(t_pv(x,y)*h_w*x_m_d*p_cross+cp_water*...
        mass_rate(x,y)*t_up)/(cp_water*mass_rate(x,y)+...
        h_w*x_m_d*p_cross);
    end
    for y = 2:y_m-1
        t_pv(x,y) = (I_m+h_a*t_a*(y_m_d*z_d+x_m_d*y_m_d)+lamda_pv*z_d*...
            ((t_pv(x,y+1)+t_pv(x,y-1))*x_m_d/y_m_d+t_pv(x-1,y)...
            *y_m_d/x_m_d)+channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
            /(h_a*(y_m_d*z_d+x_m_d*y_m_d)+lamda_pv*z_d*...
            (2*x_m_d/y_m_d+y_m_d/x_m_d)+channel_flag(x,y)*h_w*x_m_d*p_cross);
        if channel_flag(x,y)==1
            switch 3*u(x,y)+v(x,y)
                case 3
                    t_up = t_water(x-1,y);
                case -3
                    t_up = t_water(x+1,y);
                case 1
                    t_up = t_water(x,y-1);
                case -1
                    t_up = t_water(x,y+1);
                case 4
                    t_up = (t_water(x-1,y)+t_water(x,y-1))/2;
                case 2
                    t_up = (t_water(x-1,y)+t_water(x,y+1))/2;
                case -2
                    t_up = (t_water(x+1,y)+t_water(x,y))/2;
                case -4
                    t_up = (t_water(x+1,y)+t_water(x,y+1))/2;
            end
            if Re(x,y)>2300
                if t_pv(x,y)>t_water(x,y)
                    Nu= 0.023*Re(x,y)^0.8*Pr^0.4;
                else
                    Nu= 0.023*Re(x,y)^0.8*Pr^0.3;
                end
            else
                if t_pv(x,y)>t_water(x,y)
                    Nu= 5;
                else
                    Nu= 4.5;
                end
            end
            h_w = Nu*lamda_water/d_h;
            t_water(x,y)=(t_pv(x,y)*h_w*x_m_d*p_cross+cp_water*...
                mass_rate(x,y)*t_up)/(cp_water*mass_rate(x,y)+...
                h_w*x_m_d*p_cross);
        end
    end
    y=y_m;
    t_pv(x,y) = (I_m+h_a*t_a*(y_m_d*z_d+x_m_d*y_m_d)+...
        lamda_pv*z_d*(t_pv(x-1,y)*y_m_d/x_m_d+t_pv(x,y-1)*x_m_d/y_m_d)+...
        channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
        /(h_a*(y_m_d*z_d+x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+y_m_d/x_m_d)...
        +channel_flag(x,y)*h_w*x_m_d*p_cross);
    if channel_flag(x,y)==1
    switch 3*u(x,y)+v(x,y)
        case 3
            t_up = t_water(x-1,y);
        case -3
            t_up = t_water(x+1,y);
        case 1
            t_up = t_water(x,y-1);
        case -1
            t_up = t_water(x,y+1);
        case 4
            t_up = (t_water(x-1,y)+t_water(x,y-1))/2;
        case 2
            t_up = (t_water(x-1,y)+t_water(x,y+1))/2;
        case -2
            t_up = (t_water(x+1,y)+t_water(x,y))/2;
        case -4
            t_up = (t_water(x+1,y)+t_water(x,y+1))/2;
    end
    if Re(x,y)>2300
        if t_pv(x,y)>t_water(x,y)
            Nu= 0.023*Re(x,y)^0.8*Pr^0.4;
        else
            Nu= 0.023*Re(x,y)^0.8*Pr^0.3;
        end
    else
        if t_pv(x,y)>t_water(x,y)
            Nu= 5;
        else
            Nu= 4.5;
        end
    end
    h_w = Nu*lamda_water/d_h;
    t_water(x,y)=(t_pv(x,y)*h_w*x_m_d*p_cross+cp_water*...
        mass_rate(x,y)*t_up)/(cp_water*mass_rate(x,y)+...
        h_w*x_m_d*p_cross);
    end
    

end

%% Results
subplot(2,1,2)
[C,h]=contourf(t_pv(3:x_m-2,:)');
clabel(C,h)

t_avg=sum(sum(t_water))/sum(sum(channel_flag));
t_out=2*t_avg-t_water_ini
thermal_eff=cp_water*(2*t_avg-2*t_water_ini)*a_cross*1000*velo_amp/I/(numel(t_water)*x_m_d*y_m_d)