%%This script is to calculate the mismatch loss of the PV arrays
%PV temperature field
%P141
%T_pv=[t_pv t_pv t_pv t_pv t_pv t_pv];
%N
%T_pv=[t_pv fliplr(t_pv) t_pv fliplr(t_pv) t_pv fliplr(t_pv) t_pv fliplr(t_pv) t_pv fliplr(t_pv) t_pv fliplr(t_pv)];
%N12N
T_pv=t_pv;
%PV placement
T_Array=mat2cell(T_pv,[3 14*ones(1,13) 2],14*ones(1,6));
%Temperature of each module
T_module=ones(size(T_Array)-[2 0]);
for i=1:size(T_module,1)
    for j=1:size(T_module,2)
        %Average temperature of all elements
        T_module(i,j)=sum(sum(T_Array{i+1,j}))/numel(T_Array{i+1,j});
    end
end
%%calculate properties of each module
%NOCT coditions
Area_noct=2.598*0.37;
T_noct=25;
I_noct=1000;
Vop_noct=37.6;
Isc_noct=4.52*(I/I_noct)*(196*x_m_d*y_m_d/Area_noct);
Pmax_noct=115*(I/I_noct)*(196*x_m_d*y_m_d/Area_noct);
FF_noct=Pmax_noct/(Vop_noct*Isc_noct);
c_Vop=-0.28/100;
c_Isc=0.0006/100;
c_Pmax=-0.38/100;
%Initialization
Vop=ones(size(T_module));
Isc=Vop;
FF=Vop;
Pmax=Vop;
%Caculation
for i=1:size(T_module,1)
    for j=1:size(T_module,2)
        Vop(i,j)=Vop_noct*exp(c_Vop*(T_module(i,j)-T_noct));
        Isc(i,j)=Isc_noct*exp(c_Isc*(T_module(i,j)-T_noct));
        Pmax(i,j)=Pmax_noct*exp(c_Pmax*(T_module(i,j)-T_noct));
        FF(i,j)=Pmax(i,j)/(Vop(i,j)*Isc(i,j));
    end
end

%%Model and calculate mismatch loss
syms C
eqn= C^2/((1+C)*(C+log(1+C)))==mean(mean(FF));
C_overal=solve(eqn,C);
%Connection in array
Array_original=cell(size(T_module,1),size(T_module,2));
for i=1:size(Array_original,1)
    for j=1:size(Array_original,2)
        Array_original{i,j}=[i,j];
    end
end

%There are several type of patterns
pattern=1:5;
%Array is a cell storing the index for the module
for p=1:length(pattern)
    switch pattern(p)
        case 1
            %A:Vertically parallel
            Array=Array_original;
        case 2
            %B:Horizontally parallel
            Array=Array_original';
        case 3
            %C:All in sieres
            Array=cell(size(Array_original,1)*size(Array_original,2),1);
            k=1;
            for i=1:size(Array_original,1)
                for j=1:size(Array_original,2)
                    Array{k,1}=Array_original{i,j};
                    k=k+1;
                end
            end
        case 4
        %D:2nd Vertically
            m=2;
            Array=cell(size(Array_original,1)*size(Array_original,2)/m,m);
            for l=1:m
                k=1;
                for i=1:size(Array_original,1)
                    for j=size(Array_original,2)/m*(l-1)+(1:size(Array_original,2)/m)
                    Array{k,l}=Array_original{i,j};
                    k=k+1;
                    end
                end
            end
        case 5
            m=3;
            Array=cell(size(Array_original,1)*size(Array_original,2)/m,m);
            for l=1:m
                k=1;
                for i=1:size(Array_original,1)
                    for j=size(Array_original,2)/m*(l-1)+(1:size(Array_original,2)/m)
                    Array{k,l}=Array_original{i,j};
                    k=k+1;
                    end
                end
            end
    end

L=size(Array,1);
M=size(Array,2);

%Rearrange the module as array defines
Array_Vop=ones(L,M);
Array_Isc=ones(L,M);
Array_FF=ones(L,M);
Array_Pmax=ones(L,M);
for i=1:L
    for j=1:M
        Array_Vop(i,j)=Vop(Array{i,j}(1),Array{i,j}(2));
        Array_Isc(i,j)=Isc(Array{i,j}(1),Array{i,j}(2));
        Array_FF(i,j)=FF(Array{i,j}(1),Array{i,j}(2));
        Array_Pmax(i,j)=Pmax(Array{i,j}(1),Array{i,j}(2));
    end
    
end

P_ideal=sum(Array_Pmax);
omiga=P_ideal/mean(P_ideal);
P_MML=0;
for q=1:M
    syms C
    eqn= C^2/((1+C)*(C+log(1+C)))==mean(Array_FF(:,q));
    Cq=solve(eqn,C);
    sigma_I=std(Array_Isc(:,q))/mean(Array_Isc(:,q));
    P_MML=P_MML+omiga(q)*(Cq+2)/2*sigma_I^2;
end
pattern(p)
P_MML=(1-1/L)/M*P_MML+(C_overal+2)/2*std(Array_Vop(:))^2/L*[1-1/M]
clear Array
end
'P_ideal'
sum(sum(Array_Pmax))