%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Numerical modelling and design of PVT modules               %%
%                      Semester Project                                   %
%                       Xingliang Fang                                    %
%             Chair of Architecture and Building Systems                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Electric simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization

%PV temperature field from thermal simulation
if channel_type=='P141'
    T_pv=[t_pv t_pv t_pv t_pv t_pv t_pv];
elseif channel_type=='N12N'
    T_pv=t_pv;
elseif channel_type=='MEBU'
    T_pv=t_pv;
else
    T_pv=[t_pv fliplr(t_pv) t_pv fliplr(t_pv) t_pv fliplr(t_pv) t_pv fliplr(t_pv) t_pv fliplr(t_pv) t_pv fliplr(t_pv)];     
end

%% Remesh

%Dimensions of each PV cell
[x_m y_m]=size(T_pv);
if pv_type=='MONO'
    pv_cell=[14 14];
elseif pv_placement=='long'
    pv_cell=[182 14];
else
    pv_cell=[14 84];
end
T_Array=mat2cell(T_pv,[3 pv_cell(1)*ones(1,floor(x_m/pv_cell(1))) 2],pv_cell(2)*ones(1,floor(y_m/pv_cell(2))));
nb_el_cell=pv_cell(1)*pv_cell(2); %number of element per cell

%Re-calculate NOCT conditions
Vop_noct=V_op_noct;
Isc_noct=I_sc_noct*(I/I_noct)*(nb_el_cell*x_m_d*y_m_d/Area_noct);
Pmax_noct=P_max_noct*(I/I_noct)*(nb_el_cell*x_m_d*y_m_d/Area_noct);
FF_noct=Pmax_noct/(Vop_noct*Isc_noct);

%Temperature of each PV cell
T_module=ones(size(T_Array)-[2 0]);
for i=1:size(T_module,1)
    for j=1:size(T_module,2)
        %Average temperature of all elements
        T_module(i,j)=sum(sum(T_Array{i+1,j}))/numel(T_Array{i+1,j});
    end
end

%% Calculate properties of each module
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

if pv_type=='MONO'
    if pv_connect=='s'
        Array=cell(size(Array_original,1)*size(Array_original,2),1);
        k=1;
        for i=1:size(Array_original,1)
            for j=1:size(Array_original,2)
                Array{k,1}=Array_original{i,j};
                k=k+1;
            end
        end
    elseif pv_placement=='long'
        Array=Array_original;
    else
        Array=Array_original';
    end
elseif pv_placement=='long'
    if pv_connect=='s'
        Array=Array_original';
    else
        Array=Array_original;
    end
else
    if pv_connect=='s'
        Array=Array_original;
    else
        Array=Array_original';
    end
end


%Array is a cell storing the index for the module
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

if pv_type=='CIGS'
    if pv_connection=='p'
        P_ideal=Array_Pmax;
    end
else
    P_ideal=sum(Array_Pmax);
end
omiga=P_ideal/mean(P_ideal);
P_MML=zeros(1,2);
for q=1:M
    syms C
    eqn= C^2/((1+C)*(C+log(1+C)))==mean(Array_FF(:,q));
    Cq=solve(eqn,C);
    sigma_I=std(Array_Isc(:,q))/mean(Array_Isc(:,q));
    P_MML(1)=P_MML(1)+omiga(q)*(Cq+2)/2*sigma_I^2;
end
P_MML(1)=(1-1/L)/M*P_MML(1)+(C_overal+2)/2*std(Array_Vop(:))^2/L*[1-1/M];
clear Array
P_sum=sum(sum(Array_Pmax));
ele_eff=(1-P_MML(1))*P_sum/(I*x_d*y_d*6)