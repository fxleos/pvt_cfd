u=zeros(x_m,y_m);
v=zeros(x_m,y_m);
m=zeros(x_m,y_m);
switch channel_type
    case 'P'
    %Parallel
        start_p=[1,3;1,5];
        end_p=[x_m,3;x_m,5];
        m_in=[1;1];
        velo_in=[1,0;1,0];
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
