%%This script is to calculate the temperature of water

if channel_flag(x,y)==1
    switch [sign(u(x,y)),sign(v(x,y))]
        case [1,0]
            t_up = t_water(x-1,y);
        case [-1,0]
            t_up = t_water(x+1,y);
        case [0,1]
            t_up = t_water(x,y-1);
        case [0,-1]
            t_up = t_water(x,y+1);
        case [1,1]
            t_up = (t_water(x-1,y)+t_water(x,y-1))/2;
        case [1,-1]
            t_up = (t_water(x-1,y)+t_water(x,y+1))/2;
        case [-1,1]
            t_up = (t_water(x+1,y)+t_water(x,y))/2;
        case [-1,-1]
            t_up = (t_water(x+1,y)+t_water(x,y+1))/2;
    end
    t_water(x,y)=(t_pv(x,y)*h_w*x_m_d*p_cross+cp_water*...
        mass_rate(x,y)*t_up)/(cp_water*mass_rate(x,y)+...
        h_w*x_m_d*p_cross);
end