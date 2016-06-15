function [ u , v ] = velocity( velo_in,channel_flag,in_let )
% This function is to define the velocity based on geometry
%   Detailed explanation goes here
[x_m,y_m]=size(channel_flag);
u_t=zeros(x_m+2,y_m+2);
v_t=zeros(x_m+2,y_m+2);
c_f_t=zeros(x_m+2,y_m+2);
c_f_t(2:(x_m+1),2:(y_m+1))=channel_flag;
u_t(in_let(1)+1,in_let(2)+1)=velo_in(1);
v_t(in_let(1)+1,in_let(2)+1)=velo_in(2);
x=in_let(1)+1;
y=in_let(2)+1;

while x<x_m+2 && y<y_m+2
    switch 2*u_t(x,y)+v_t(x,y)
        case 1
            switch 4*c_f_t(x+1,y)+2*c_f_t(x,y+1)+c_f_t(x-1,y)
                case 1
                    u_t(x-1,y)=u_t(x-1,y)-1;
                    x=x-1;
                case 2
                    v_t(x,y+1)=v_t(x,y+1)+1;
                    y=y+1;
                case 4
                    u_t(x+1,y)=u_t(x+1,y)+1;
                    x=x+1;
                case 0
                    x=x_m+2;
            end
        case 2
            switch 4*c_f_t(x,y-1)+2*c_f_t(x+1,y)+c_f_t(x,y+1)
                case 1
                    v_t(x,y+1)=v_t(x,y+1)+1;
                    y=y+1;
                case 2
                    u_t(x+1,y)=u_t(x+1,y)+1;
                    x=x+1;
                case 4
                    v_t(x,y-1)=v_t(x,y-1)-1;
                    y=y-1;
                case 0
                    x=x_m+2;
            end
        case -1
            switch 4*c_f_t(x-1,y)+2*c_f_t(x,y-1)+c_f_t(x+1,y)
                case 1
                    u_t(x+1,y)=u_t(x+1,y)+1;
                    x=x+1;
                case 2
                    v_t(x,y-1)=v_t(x,y-1)-1;
                    y=y-1;
                case 4
                    u_t(x-1,y)=u_t(x-1,y)-1;
                    x=x-1;
                case 0
                    x=x_m+2;
            end
        case -2
            switch 4*c_f_t(x,y-1)+2*c_f_t(x-1,y)+c_f_t(x,y+1)
                case 1
                    v_t(x,y+1)=v_t(x,y+1)+1;
                    y=y+1;
                case 2
                    u_t(x-1,y)=u_t(x-1,y)-1;
                    x=x-1;
                case 4
                    v_t(x,y-1)=v_t(x,y-1)-1;
                    y=y-1;
                case 0
                    x=x_m+2;
            end
            
    end

    
end

u=u_t(2:x_m+1,2:y_m+1);
v=v_t(2:x_m+1,2:y_m+1);
end


