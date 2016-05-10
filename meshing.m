%%This script is to draw geometry and do meshing
% Meshing Demisions
if channel_type=='P141'
    y_m=14;      %Element number in x direction
else 
    y_m=7;      %Element number in x direction
end
y_m_d=y_d/y_m;  %Element length in x direction
x_m=round(x_d/y_m_d);        %Element number in y direction with similar element length
x_m_d=x_d/x_m; %Element length in y direction

channel_flag=zeros(x_m,y_m);

switch channel_type
    case 'P'
    %Parallel
    channel_flag(:,3)=1;
    channel_flag(:,5)=1;
    case 'P141'
        %Parallel 1-4-1
        channel_flag(1:6,7)=1;
        channel_flag((x_m-5):x_m,8)=1;
        channel_flag([7,x_m-6],6:9)=1;
        channel_flag([8,9,x_m-7,x_m-8],[6,9])=1;
        channel_flag([9,x_m-8],3:6)=1;
        channel_flag([9,x_m-8],9:12)=1;
        channel_flag(10:(x_m-9),[3,5,10,12])=1;
    case 'U'
    %U
    channel_flag(1:7,[1,7])=1;
    channel_flag(7,[2,6])=1;
    channel_flag(7:x_m-6,[3,5])=1;
    channel_flag(x_m-6,4)=1;
    case 'N'
    %N
    channel_flag(1:7,1)=1;
    channel_flag(7:(x_m-6),[2,4,6])=1;
    channel_flag((x_m-6):x_m,y_m)=1;
    channel_flag(x_m-6,3)=1;
    channel_flag(7,5)=1;
end
