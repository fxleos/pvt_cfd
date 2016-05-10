function [ u , v , m] = velocity_water( velo_in,m_in,channel_flag,start_p,end_p )
% This function is to define the velocity based on geometry
%   Detailed explanation goes here
[x_m,y_m]=size(channel_flag);
u_t=zeros(x_m+2,y_m+2);
v_t=zeros(x_m+2,y_m+2);
m_t=zeros(x_m+2,y_m+2);
c_f_t=zeros(x_m+2,y_m+2);
c_f_t(2:(x_m+1),2:(y_m+1))=channel_flag;
u_t(start_p(1)+1,start_p(2)+1)=velo_in(1);
v_t(start_p(1)+1,start_p(2)+1)=velo_in(2);
m_t(start_p(1)+1,start_p(2)+1)=m_in;
x=start_p(1)+1;
y=start_p(2)+1;

while ~(x==end_p(1)+1 && y==end_p(2)+1)
    switch 2*u_t(x,y)+v_t(x,y)
        case 1
            switch 4*c_f_t(x+1,y)+2*c_f_t(x,y+1)+c_f_t(x-1,y)
                case 1
                    u_t(x-1,y)=u_t(x-1,y)-1;
                    m_t(x-1,y)=m_t(x-1,y)+m_in;
                    x=x-1;
                case 2
                    v_t(x,y+1)=v_t(x,y+1)+1;
                    m_t(x,y+1)=m_t(x,y+1)+m_in;
                    y=y+1;
                case 4
                    u_t(x+1,y)=u_t(x+1,y)+1;
                    m_t(x+1,y)=m_t(x+1,y)+m_in;
                    x=x+1;
                case 0
                    x=end_p(1)+1; y=end_p(2)+1;
            end
        case 2
            switch 4*c_f_t(x,y-1)+2*c_f_t(x+1,y)+c_f_t(x,y+1)
                case 1
                    v_t(x,y+1)=v_t(x,y+1)+1;
                    m_t(x,y+1)=m_t(x,y+1)+m_in;
                    y=y+1;
                case 2
                    u_t(x+1,y)=u_t(x+1,y)+1;
                    m_t(x+1,y)=m_t(x+1,y)+m_in;
                    x=x+1;
                case 4
                    v_t(x,y-1)=v_t(x,y-1)-1;
                    m_t(x,y-1)=m_t(x,y-1)+m_in;
                    y=y-1;
                case 0
                    x=end_p(1)+1; y=end_p(2)+1;
            end
        case -1
            switch 4*c_f_t(x-1,y)+2*c_f_t(x,y-1)+c_f_t(x+1,y)
                case 1
                    u_t(x+1,y)=u_t(x+1,y)+1;
                    m_t(x+1,y)=m_t(x+1,y)+m_in;
                    x=x+1;
                case 2
                    v_t(x,y-1)=v_t(x,y-1)-1;
                    m_t(x,y-1)=m_t(x,y-1)+m_in;
                    y=y-1;
                case 4
                    u_t(x-1,y)=u_t(x-1,y)-1;
                    m_t(x-1,y)=m_t(x-1,y)+m_in;
                    x=x-1;
                case 0
                    x=end_p(1)+1; y=end_p(2)+1;
            end
        case -2
            switch 4*c_f_t(x,y-1)+2*c_f_t(x-1,y)+c_f_t(x,y+1)
                case 1
                    v_t(x,y+1)=v_t(x,y+1)+1;
                    m_t(x,y+1)=m_t(x,y+1)+m_in;
                    y=y+1;
                case 2
                    u_t(x-1,y)=u_t(x-1,y)-1;
                    m_t(x-1,y)=m_t(x-1,y)+m_in;
                    x=x-1;
                case 4
                    v_t(x,y-1)=v_t(x,y-1)-1;
                    m_t(x,y-1)=m_t(x,y-1)+m_in;
                    y=y-1;
                case 0
                    x=end_p(1)+1;y=end_p(2)+1;
            end
            
    end

    
end

u=u_t(2:x_m+1,2:y_m+1);
v=v_t(2:x_m+1,2:y_m+1);
m=m_t(2:x_m+1,2:y_m+1);
end


