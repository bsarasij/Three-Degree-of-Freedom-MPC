clear
close all
clc
load inp_out_db_26_5.mat
load modparms_18_5_ms.mat

modparms = modparms_18_5_ms;
modparms.kern = 'tric';

ctrl_data = inp_out_db_26_5;
ctrl_data_ms = ctrl_data-mean(ctrl_data);
narx_reg_struct = [4 4 4 4 6 1 1];
Tsamp = 1;
[Ydb,phi_db]=modmkreg(ctrl_data_ms,narx_reg_struct,Tsamp); % Create Regressor for each day
KS=[650 750];

%%Generating the i varying regressor vector for MoD-estimation
n_out=size(narx_reg_struct,1);
n_in=(size(narx_reg_struct,2)-n_out)/2;

NA=narx_reg_struct(:,1:n_out); % for our case,  n_out=2
NB=narx_reg_struct(:,n_out+1:n_out+n_in); % for our case, ; n_in=3
NK=narx_reg_struct(:,n_out+n_in+1:n_out+2*n_in);
ll=1;

na = NA(ll,:); nb = NB(ll,:);
regs = sum(na)+sum(nb); % reg=2
        
Y_mod = modparms.Z.Y(ll,:);
X_mod = modparms.Z.X(1:regs,:,ll);
Scale=modparms.Z.M;
GOF = {modparms.gof,modparms.pen};
guimm=modparms.minmeth;

SS = Scale(1:regs,1:regs,ll);


figure;
for iter = 1:size(phi_db,2)
    mod_regvec = phi_db(:,iter);
    beta=zeros(sum(NA(1,:))+sum(NB(1,:))+1,n_out);
    var=[];
    ind2=[];
        x_nn = mod_regvec;
        if n_out > 1
            ind2 = 1:regs;
        end

        [beta_ind,~,Si] = locpol(Y_mod,X_mod,x_nn,1,SS,modparms.kern,GOF,guimm,1./var,KS,ind2);
        beta(:,ll)=beta_ind;
        if isempty(modparms.Z.var)
            var = Si.var;
        end
   
    beta1=beta(2:end,:);
    [id_A,id_B]=beta2idpoly(beta1,n_out,n_in,NA,NB);
    idpol=idpoly(id_A,id_B,'Ts',Tsamp);
    ss_mod_orig=ss(idpol);
    ss_mod_orig.InputDelay=NK'-1;
    ss_model=ssDel2ABCD(ss_mod_orig);
    step_response = step(ss_model);
    subplot(1,3,1);plot(step_response(:,1,1));hold on;xlim([0 500]);subplot(1,3,2);plot(step_response(:,1,2));hold on;xlim([0 500]);subplot(1,3,3);plot(step_response(:,1,3));hold on;xlim([0 500]);
    disp(iter);
    pause(0.1);
end
