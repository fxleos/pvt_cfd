 if channel_flag(x,y)==1
                if velo(x,y)>0
                    if channel_flag(x+1,y)==1
                        t_d1=t_water(x+1,y);
                        if x<x_m-1
                            t_d2=t_water(x+2,y);
                        else
                            t_d2=t_d1+0.01;
                        end
                    elseif channel_flag(x,y+1)==1
                        t_d1=t_water(x,y+1);
                        if x<x_m-1
                            t_d2=t_water(x,y+2);
                        else
                            t_d2=t_d1+0.01;
                        end
                    end
                    if channel_flag(x-1,y)==1
                        t_u1=t_water(x-1,y);
                        if x>2
                            t_u2=t_water(x-2,y);
                        else
                            t_u2=t_u1-0.01;
                        end
                    elseif channel_flag(x,y-1)==1
                        t_u1=t_water(x,y-1);
                        if x>2
                            t_u2=t_water(x,y-2);
                        else
                            t_u2=t_u1-0.01;
                        end
                    end
                elseif velo(x,y)<0
                    if channel_flag(x-1,y)==1
                        t_d1=t_water(x-1,y);
                        if x>2
                            t_d2=t_water(x-2,y);
                        else
                            t_d2=t_d1+0.01;
                        end
                    elseif channel_flag(x,y-1)==1
                        t_d1=t_water(x,y-1);
                        if x>2
                            t_d2=t_water(x,y-2);
                        else
                            t_d2=t_d1+0.01;
                        end
                    end
                    if channel_flag(x+1,y)==1
                        t_u1=t_water(x+1,y);
                        if x<x_m-1
                            t_u2=t_water(x+2,y);
                        else
                            t_u2=t_u1-0.01;
                        end
                    elseif channel_flag(x,y+1)==1
                        t_u1=t_water(x,y+1);
                        if x<x_m-1
                            t_u2=t_water(x,y+2);
                        else
                            t_u2=t_u1-0.01;
                        end
                    end
                end
                t_water(x,y)=t_pv(x,y)-cp_water*mass_rate(x,y)*(t_d1+t_d2-t_u1-t_u2)/...
                    (2*h_w*x_m_d*p_cross); %the p_cross*x_m_d is not accurate
                        
            end