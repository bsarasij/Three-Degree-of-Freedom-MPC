function uk=control_gen(pH_k,CO_2_k_1,Dm_curr, filt_ref, time_i, Unmeas_DistModel,degree,slack_var,anticipation, Weights,limits,setpoints,p,m,nud,fa,D_forecast_filt,D_forecast_unfilt,Dm_filt_minus1,Dm_filt_minus2,mld_copy,model_choice,ts)

persistent modparms narx_reg_struct big_reg_pH big_reg_CO_2 big_reg_Temp big_reg_Rad x_hat_k_given_k x_hat_k_given_k_filt 


% Load modparams
if isempty(modparms)
    out_modparms=load('modparms_non_ms.mat');
    modparms=out_modparms.modparms;
end

% if ~exist(modparms, 'var')
%     load('', 'modparms');
% end

modparms.kern='tric';
if isempty(narx_reg_struct)
    out_narx=load('narx_reg_struct.mat');
    narx_reg_struct=out_narx.narx_reg_struct;
end



%%Generating the i varying regressor vector for MoD-estimation
n_out=size(narx_reg_struct,1);
n_in=(size(narx_reg_struct,2)-n_out)/2;

NA=narx_reg_struct(:,1:n_out); % for our case,  n_out=2
NB=narx_reg_struct(:,n_out+1:n_out+n_in); % for our case, ; n_in=3
NK=narx_reg_struct(:,n_out+n_in+1:n_out+2*n_in);

if time_i<=1
    big_reg_pH = zeros(NA(1,1)+1,1); % Stores [y_k y_k-1 y_k-2... y_k-na]
    big_reg_CO_2=zeros(NB(1,1)+NK(1,1)-1,1); %[u_k-1 u_k-2 u_k-3... u_k-nb-nk+1]
    big_reg_Temp=zeros(NB(1,2)+NK(1,2),1); %[Dm_k Dm_k-1 Dm_k-2... Dm_k-nb-nk]
    big_reg_Rad=zeros(NB(1,3)+NK(1,3),1);
    x_init=zeros(sum(NA(1,:))+sum(NB(1,:))+2*n_out-1,1); % Defining Intitial (Steady-State) Conditions
    x_hat_k_given_k=x_init;
    x_hat_k_given_k_filt=x_init;
end
Temp_k=Dm_curr(1);
Rad_k=Dm_curr(2);
big_reg_pH=circshift(big_reg_pH,1);
big_reg_pH(1,:)=pH_k;

big_reg_CO_2=circshift(big_reg_CO_2,1);
big_reg_CO_2(1,:)=CO_2_k_1;

big_reg_Temp=circshift(big_reg_Temp,1);
big_reg_Temp(1,:)=Temp_k;

big_reg_Rad=circshift(big_reg_Rad,1);
big_reg_Rad(1,:)=Rad_k;

%%Assigning the values needed for control calculation
x_hat_k_minus_1_given_k_minus_1 = x_hat_k_given_k;
x_hat_k_minus_1_given_k_minus_1_filt=x_hat_k_given_k_filt;
mod_regvec = [big_reg_pH(2:end);big_reg_CO_2(NK(1,1):NB(1,1)+NK(1,1)-1);big_reg_Temp(1+NK(1,2):NB(1,2)+NK(1,2));big_reg_Rad(1+NK(1,3):NB(1,3)+NK(1,3))]; % Regressor vector for MoD
y_plant_measurement=big_reg_pH(1);
uk_minus_1=big_reg_CO_2(1);
uk_minus_2=big_reg_CO_2(2);
delta_uk_minus_1=uk_minus_1-uk_minus_2;
Dm_minus1=[big_reg_Temp(2);big_reg_Rad(2)];
Dm_minus2=[big_reg_Temp(3);big_reg_Rad(3)];
deltak_minus_1=0;%Set to zero as currently no categorical variable
zk_minus_1=0;
delta_deltak_minus_1=0;
delta_zk_minus_1=0;

%% Generating Predictive Models -MoD/ARX

if model_choice==1
    KS=[650 750];
    beta=zeros(sum(NA(1,:))+sum(NB(1,:))+1,n_out);
    var=[];
    ind2=[];
    for ll=1:n_out  % for our case, n_out=2 (# of output)
        na = NA(ll,:); nb = NB(ll,:);
        regs = sum(na)+sum(nb); % reg=2
        x_nn = mod_regvec;
        if n_out > 1
            ind2 = 1:regs;
        end

        Y_mod = modparms.Z.Y(ll,:);
        X_mod = modparms.Z.X(1:regs,:,ll);
        Scale=modparms.Z.M;
        GOF = {modparms.gof,modparms.pen};
        guimm=modparms.minmeth;
        
        SS = Scale(1:regs,1:regs,ll);
        [beta_ind,~,Si] = locpol(Y_mod,X_mod,x_nn,1,SS,modparms.kern,GOF,guimm,1./var,KS,ind2);
        beta(:,ll)=beta_ind;
        if isempty(modparms.Z.var)
            var = Si.var;
        end
    end
    beta1=beta(2:end,:);
    [id_A,id_B]=beta2idpoly(beta1,n_out,n_in,NA,NB);
    idpol=idpoly(id_A,id_B,'Ts',ts);
    ss_mod_orig=ss(idpol);
    ss_mod_orig.InputDelay=NK'-1;
    ss_model=ssDel2ABCD(ss_mod_orig);

elseif model_choice==0
    arx_model=modparms.arxmodel;
    id_A_arx=arx_model.A;
    id_B_arx=arx_model.B;
    idpol=idpoly(id_A_arx,id_B_arx,'Ts',ts);
    ss_mod_orig=ss(idpol);
    ss_mod_orig.InputDelay=NK'-1;
    ss_model=ssDel2ABCD(ss_mod_orig);
end
%%
Ac=ss_model.A;
Cc=ss_model.C;
Bc=ss_model.B(:,1);
Bd_orig=ss_model.B(:,2:end);
Bd=Bd_orig;

A1=Ac;
B1=Bc;
Cmod=Cc;
B2=zeros(size(A1,1),size(mld_copy.E2,2));
B3=zeros(size(A1,1),size(mld_copy.E3,2));

sv=step(ss_model);
subplot(1,3,1);plot(sv(:,:,1));hold on;subplot(1,3,2);plot(sv(:,:,2));hold on;subplot(1,3,3);plot(sv(:,:,3));hold on;drawnow;

MldModel=struct('A',A1,'B1',B1,'B2',B2,'B3',B3,'Bd',Bd,'C',Cmod,'E1',mld_copy.E1,'E2',mld_copy.E2,'E3',mld_copy.E3,'E4',mld_copy.E4,'E5',mld_copy.E5,'Ed',mld_copy.Ed);
[HybMpcCon] = HybMpcCon_code(MldModel,Unmeas_DistModel,Weights,limits,setpoints,p,m,nud,fa,slack_var,anticipation);


[uk,x_hat_k_given_k,x_hat_k_given_k_filt]= Hyb_Mpc44_bioreactor(x_hat_k_minus_1_given_k_minus_1,y_plant_measurement,uk_minus_1,delta_uk_minus_1,Dm_minus1,Dm_minus2,Dm_filt_minus1,Dm_filt_minus2,D_forecast_filt,D_forecast_unfilt,filt_ref,deltak_minus_1,zk_minus_1,delta_deltak_minus_1,delta_zk_minus_1,x_hat_k_minus_1_given_k_minus_1_filt, HybMpcCon, degree,slack_var);

end