function [HybMpcCon]=HybMpcCon_code(MldModel,DistModel,Weights,limits,setpoints,p,m,nud,fa,slack_var,anticipation)
Fa = diag(fa);                 % See Equation 41                       
Fb=(diag(fa)^2)/(eye(size(DistModel.Aw))+DistModel.Aw-DistModel.Aw*diag(fa));  % See Equation 42-43
Kf = [zeros(size(MldModel.A,1),size(MldModel.C,1));Fb; Fa];   % See Equation 40


A_aug = [MldModel.A,  zeros(size(MldModel.A,1),size(DistModel.Aw,2)), zeros(size(MldModel.A,1),size(MldModel.C,1));...   
        zeros(size(DistModel.Aw,1),size(MldModel.A,2)), DistModel.Aw, zeros(size(DistModel.Aw,1),size(MldModel.C,1));...
        MldModel.C*MldModel.A,  DistModel.Aw,  eye(size(MldModel.C,1))];
B1_aug = [MldModel.B1; zeros(size(DistModel.Aw,1),size(MldModel.B1,2));  MldModel.C*MldModel.B1];
B2_aug = [MldModel.B2; zeros(size(DistModel.Aw,1), size(MldModel.B2,2)); MldModel.C*MldModel.B2];
B3_aug = [MldModel.B3; zeros(size(DistModel.Aw,1),  size(MldModel.B3,2)); MldModel.C*MldModel.B3];
Bd_aug = [MldModel.Bd;   zeros(size(DistModel.Aw,1),size(MldModel.Bd,2)); MldModel.C*MldModel.Bd]; 
Bw_aug = [zeros(size(MldModel.A,1),size(DistModel.Aw,2)); eye(size(DistModel.Aw,1)); eye(size(MldModel.C,1))];  
C_aug = [zeros(size(MldModel.C,1),size(MldModel.A,2)) zeros(size(MldModel.C,1),size(DistModel.Aw,2)) eye(size(MldModel.C,1))];

[~,nu_aug]=size(B1_aug);
[~,nd_aug]=size(B2_aug);
[~,nz_aug]=size(B3_aug);
% [ny,~] = size(C_aug);
%==================weight matrices for u,y& z==============================
dig_wy=repmat(Weights.wy,p,1);
Wy=diag(dig_wy);                            % weight on output error

dig_wu=repmat(Weights.wu,m,1);
Wu=diag(dig_wu);                            % weight on imove suppression

dig_wdu=repmat(Weights.wdu,m,1);
Wdu=diag(dig_wdu);

dig_wd=repmat(Weights.wd,p,1);
Wd=diag(dig_wd);

dig_wz=repmat(Weights.wz,p,1);
Wz=diag(dig_wz);


diag_wslack = repmat(Weights.wslack,p,1);
Q_slack = diag(diag_wslack);
%===================================generates H matrix=====================

[Qhat,Sb,phi,H1,H2,H3,Hd,H11,H21,H31,Hd1,Epsilon_1,Epsilon_2,Epsilon_3,Epsilon_4,Epsilon_5,Epsilon_d,Epsilon_41,Epsilon_42,Epsilon_43,Epsilon_4d,Ru,Ru0]= MLD_matrices(p,m,A_aug,B1_aug,B2_aug,B3_aug,Bd_aug,C_aug,MldModel,Wu,Wz,Wy,Wdu,Wd,Q_slack,slack_var,anticipation);
H=2*Qhat;


Yu=repmat(limits.yu,p,1);
Yl=repmat(limits.yl,p,1);
Uu=repmat(limits.uu,m,1);
Ul=repmat(limits.ul,m,1);
delta_Uu=repmat(limits.deluu,m,1);
delta_Ul=repmat(limits.delul,m,1);

Delta_up=repmat(limits.delta_up,p,1);
Delta_lo=repmat(limits.delta_lo,p,1);
Z_up=repmat(limits.z_up,p,1);
Z_lo=repmat(limits.z_lo,p,1);
slack_lo = repmat(limits.slack_min,p,1);
slack_up = repmat(limits.slack_max,p,1);
 
if slack_var ==0
lb=[Ul;Delta_lo;Z_lo];
ub=[Uu;Delta_up;Z_up];
elseif slack_var == 1
lb=[Ul;Delta_lo;Z_lo;slack_lo];
ub=[Uu;Delta_up;Z_up;slack_up];
end
% =========================================================================
U_sp=repmat(setpoints.u,m,1);     %U_sp1=[u_sp1 u_sp1]'
Del_sp=repmat(setpoints.del,p,1); %Del_sp1=[del_sp1 del_sp1]'
Z_sp=repmat(setpoints.z,p,1);     %Z_sp1=[z_sp1 z_sp1 z_sp1 z_sp1 z_sp1]';



nuc_aug=nu_aug-nud;
IntVars=zeros(m*nu_aug+p*nd_aug+p*nz_aug,1);
for i=1:m
    IntVars((i-1)*nu_aug+nuc_aug+1:i*nu_aug,1)=ones(nud,1); % Intiger variables
end
IntVars(m*nu_aug+1:m*nu_aug+p*nd_aug)=ones(p*nd_aug,1);
kk=1;
ll=1;
for i=1:length(IntVars)
    if IntVars(i)==1
        binary_var_index(kk)=i;
        kk=kk+1;
    else
        contivar(ll)=i;
        ll=ll+1;
    end
end
vartype=binary_var_index;

HybMpcCon=struct('A',A_aug,'B1',B1_aug,'B2',B2_aug,'B3',B3_aug,'Bd',Bd_aug,'Bw',Bw_aug,'C',C_aug,'Kf',Kf,'Yu',Yu,'Yl',Yl,'Uu',Uu,'Ul',Ul,'delta_Uu',delta_Uu,'delta_Ul',delta_Ul,...
    'phi',phi,'H1',H1,'H2',H2,'H3',H3,'Hd',Hd,'H11',H11,'H21',H21,'H31',H31,'Hd1',Hd1,'Epsilon_1', Epsilon_1, 'Epsilon_2',Epsilon_2,'Epsilon_3',Epsilon_3,'Epsilon_4',Epsilon_4,'Epsilon_5',Epsilon_5,'Epsilon_d', Epsilon_d, 'Epsilon_41',Epsilon_41,'Epsilon_42',Epsilon_42,'Epsilon_43',Epsilon_43,'Epsilon_4d',Epsilon_4d,'Ru',Ru,'Ru0',Ru0,'Wy',Wy,'Wu',Wu,'Wdu',Wdu,'Wd',Wd,'Wz',Wz,...
    'H',H,'Sb',Sb,'lb',lb,'ub',ub,'vartype',vartype,'Z_up',Z_up,'U_sp',U_sp,'Del_sp',Del_sp, 'Delta_up',Delta_up, 'Z_sp',Z_sp,'p',p,'m',m,'Q_slack',Q_slack);
% MeasDist=struct('filt',MeasDist_filt);







