%% Simulation

for i = 1:it
    x=1;y=1;
    t_pv(x,y) = (I_m+h_a*t_a*(y_m_d*z_d+2*x_m_d*y_m_d)+...
        lamda_pv*z_d*(t_pv(x,y+1)*x_m_d/y_m_d+t_pv(x+1,1)*y_m_d/x_m_d)+...
        channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
        /(h_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+y_m_d/x_m_d)...
        +channel_flag(x,y)*h_w*x_m_d*p_cross);
    %t_pv(1,1) = (I_m+h_a*t_a*(y_m_d*z_d+x_m_d*y_m_d)+lamda_pv*z_d*(t_pv(1,2)*...
    %    x_m_d/y_m_d+t_pv(2,1)*y_m_d/x_m_d))/(h_a*(y_m_d*z_d+x_m_d*y_m_d)...
    %    +lamda_pv*z_d*(x_m_d/y_m_d+y_m_d/x_m_d));
    temperture_water;
    
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
        
    temperture_water;
    end
    
    x=1;y=y_m;
    t_pv(x,y) = (I_m+h_a*t_a*(y_m_d*z_d+2*x_m_d*y_m_d)+...
        lamda_pv*z_d*(t_pv(x,y-1)*x_m_d/y_m_d+t_pv(x+1,y)*y_m_d/x_m_d)+...
        channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
        /(h_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+y_m_d/x_m_d)...
        +channel_flag(x,y)*h_w*x_m_d*p_cross);
    %t_pv(1, y_m) = (I_m+h_a*t_a*(y_m_d*z_d+x_m_d*y_m_d)+lamda_pv*z_d*(t_pv(1,y_m-1)*...
    %    x_m_d/y_m_d+t_pv(2,y_m)*y_m_d/x_m_d))/(h_a*(y_m_d*z_d+x_m_d*y_m_d)...
    %    +lamda_pv*z_d*(x_m_d/y_m_d+y_m_d/x_m_d));
    temperture_water;
    
    for x = 2:x_m-1
        y=1;
        t_pv(x,y) = (I_m+h_a*t_a*(2*x_m_d*y_m_d)+lamda_pv*z_d*((t_pv(x-1,1)...
            +t_pv(x+1,1))*y_m_d/x_m_d+t_pv(x,y+1)*x_m_d/y_m_d)+...
            channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
            /(h_a*(2*x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+2*y_m_d/x_m_d)...
            +channel_flag(x,y)*h_w*x_m_d*p_cross);
        temperture_water;
        y=y_m;
        t_pv(x,y) = (I_m+h_a*t_a*(2*x_m_d*y_m_d)+lamda_pv*z_d*((t_pv(x-1,1)...
            +t_pv(x+1,1))*y_m_d/x_m_d+t_pv(x,y-1)*x_m_d/y_m_d)+...
            channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
            /(h_a*(2*x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+2*y_m_d/x_m_d)...
            +channel_flag(x,y)*h_w*x_m_d*p_cross);
        %t_pv(x,1) = (I_m+h_a*t_a*(x_m_d*y_m_d)+lamda_pv*z_d*((t_pv(x-1,1)...
        %    +t_pv(x+1,1))*y_m_d/x_m_d+t_pv(x,2)*x_m_d/y_m_d))...
        %    /(h_a*(x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+2*y_m_d/x_m_d));
        %t_pv(x,y_m) = (I_m+h_a*t_a*(x_m_d*y_m_d)+lamda_pv*z_d*((t_pv(x-1,y_m)...
        %    +t_pv(x+1,y_m))*y_m_d/x_m_d+t_pv(x,y_m-1)*x_m_d/y_m_d))...
        %    /(h_a*(x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+2*y_m_d/x_m_d));
        temperture_water;
    end
    
    x=x_m;y=1;
    t_pv(x,y) = (I_m+h_a*t_a*(y_m_d*z_d+2*x_m_d*y_m_d)+...
        lamda_pv*z_d*(t_pv(x-1,y)*y_m_d/x_m_d+t_pv(x,y+1)*x_m_d/y_m_d)+...
        channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
        /(h_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+2*y_m_d/x_m_d)...
        +channel_flag(x,y)*h_w*x_m_d*p_cross);
    temperture_water;
    x=x_m;y=y_m;
    t_pv(x,y) = (I_m+h_a*t_a*(y_m_d*z_d+2*x_m_d*y_m_d)+...
        lamda_pv*z_d*(t_pv(x-1,y)*y_m_d/x_m_d+t_pv(x,y-1)*x_m_d/y_m_d)+...
        channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
        /(h_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*(x_m_d/y_m_d+2*y_m_d/x_m_d)...
        +channel_flag(x,y)*h_w*x_m_d*p_cross);
    temperture_water;
    
    x=x_m;
    for y = 2:y_m-1
        t_pv(x,y) = (I_m+h_a*t_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*...
            ((t_pv(x,y+1)+t_pv(x,y-1))*x_m_d/y_m_d+t_pv(x-1,y)...
            *y_m_d/x_m_d)+channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
            /(h_a*(y_m_d*z_d+2*x_m_d*y_m_d)+lamda_pv*z_d*...
            (2*x_m_d/y_m_d+y_m_d/x_m_d)+channel_flag(x,y)*h_w*x_m_d*p_cross);
        temperture_water;
    end
    
    for x=2:x_m-1
        for y =2:y_m-1
            t_pv(x,y) = (I_m+h_a*t_a*2*x_m_d*y_m_d+lamda_pv*z_d*...
                ((t_pv(x,y-1)+t_pv(x,y+1))*x_m_d/y_m_d+...
                ((t_pv(x-1,y)+t_pv(x+1,y))*y_m_d/x_m_d))+...
                channel_flag(x,y)*h_w*t_water(x,y)*x_m_d*p_cross)...
                /(h_a*2*x_m_d*y_m_d+lamda_pv*z_d*2*(x_m_d/y_m_d+y_m_d/x_m_d)...
                +channel_flag(x,y)*h_w*x_m_d*p_cross);
            temperture_water;
        
        end 
    end
    
    
end