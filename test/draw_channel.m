switch channel_type
    case 'P'
    %Parallel
    channel_flag(:,3)=1;
    channel_flag(:,5)=1;
    case 'P141'
        %Parallel 1-4-1
        channel_flag(1:7,1)=1;
        channel_flag((x_m-6):x_m,1)=1;
        channel_flag(7:9,2)=1;
        channel_flag((x_m-8):(x_m-6),2)=1;
        channel_flag([9,x_m-8],2:5)=1;
        channel_flag(10:(x_m-9),[3,5])=1;
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
