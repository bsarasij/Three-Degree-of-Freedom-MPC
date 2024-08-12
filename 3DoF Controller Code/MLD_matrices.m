
function [Qhat,Sb,phi,H1,H2,H3,Hd,H11,H21,H31,Hd1,Epsilon_1,Epsilon_2,Epsilon_3,Epsilon_4,Epsilon_5,Epsilon_d,Epsilon_41,Epsilon_42,Epsilon_43,Epsilon_4d,Ru,Ru0] = MLD_matrices(p,m,A_aug,B1_aug,B2_aug,B3_aug,Bd_aug,C_aug,MldModel,Wu,Wz,Wy,Wdu,Wd,Q_slack, slack_var,anticipation)
    
%% ===============================Dimensions Declaration =======================================================================
[~,nu_aug]=size(B1_aug);
[~,nd_aug]=size(B2_aug);
[~,nz_aug]=size(B3_aug);
[ny_aug,~]=size(C_aug);
[nx_aug, ndi_aug]=size(Bd_aug);
ne=size(MldModel.E1,1);


%% ============================================================Generation of Ru and Ru0 for Equation (21)============================================================

%=====================Ru and Ru0 matrices=========================
Ru_dummy=zeros(m*nu_aug,nu_aug);
Ru_dummy(1:nu_aug,:)= eye(nu_aug);
Ru_dummy(nu_aug+1:2*nu_aug,:) = -eye(nu_aug);
Ru=zeros(m*nu_aug,m*nu_aug);
for i=1:m
%=========================Ru===============================
    Ru((i-1)*nu_aug+1:m*nu_aug,(i-1)*nu_aug+1:i*nu_aug)=Ru_dummy(1:(m+1-i)*nu_aug,1:nu_aug);
end
%=========================Ru0==============================
Ru0=zeros(m*nu_aug,nu_aug);
Ru0(1:nu_aug,:)= eye(nu_aug);

%% ===============================phi matrix in Equation (A.3) =======================================================================
phi = zeros(p*ny_aug,nx_aug);
CApoweri = C_aug;
for i=1:p
    CApoweri=CApoweri*A_aug;
    phi((i-1)*ny_aug+1:i*ny_aug,1:nx_aug)=CApoweri;
end

%% ============================================== Generation of H1, H2, H3, Hd, H11, H21, H31, Hd1 in Equations (A.1-A.2, A.4)=========================================================

%========================H1,H2, H3, and Hd Dummy matrices==================
H1_dummy = zeros(p*ny_aug,nu_aug); %H1
H2_dummy = zeros(p*ny_aug,nd_aug); %H2
H3_dummy = zeros(p*ny_aug,nz_aug); %H3
Hd_dummy = zeros(p*ny_aug,ndi_aug); %Hd
Apoweriminus1=eye(size(A_aug));
Apoweriminus2=eye(size(A_aug));


H1 = zeros(p*ny_aug,m*nu_aug);%Hd
H2 = zeros(p*ny_aug,p*nd_aug); %H2
H3=zeros(p*ny_aug,p*nz_aug); %H3
Hd = zeros(p*ny_aug,p*ndi_aug);%Hd
H11 = zeros(p*ny_aug,nu_aug); %H11
H21= zeros(p*ny_aug,nd_aug); %H21
H31 = zeros(p*ny_aug,nz_aug); %H31
Hd1 = zeros(p*ny_aug,ndi_aug); %Hd1

for i=1:p
    
    if i==1
        CApoweriminus1B1 = C_aug*B1_aug; %H1
        CApoweriminus1B11 = C_aug*B1_aug;
        if nd_aug~=0
            CApoweriminus1B2 = C_aug*B2_aug; %H2
            CApoweriminus1B21 = C_aug*B2_aug;
        end
        if nz_aug~=0
            CApoweriminus1B3 = C_aug*B3_aug; %H3
            CApoweriminus1B31 = C_aug*B3_aug;
        end
        if ndi_aug~=0
            CApoweriminus1Bd = C_aug*Bd_aug; %Hd
            CApoweriminus1Bd1 = C_aug*Bd_aug;
        end
    else
        Apoweriminus1 = Apoweriminus1*A_aug;
        if i>2
            Apoweriminus2 = Apoweriminus2*A_aug;
        end
%===================For H1===================
        CApoweriminus1B1 = C_aug*Apoweriminus1*B1_aug-C_aug*Apoweriminus2*B1_aug;
        CApoweriminus1B11 = C_aug*Apoweriminus1*B1_aug;
      
 %===================For H2===================
  if nd_aug~=0           
            CApoweriminus1B2 = C_aug*Apoweriminus1*B2_aug-C_aug*Apoweriminus2*B2_aug;
            CApoweriminus1B21 = C_aug*Apoweriminus1*B2_aug;
  end
       
%===================For H3===================
 if nz_aug~=0         
             CApoweriminus1B3 = C_aug*Apoweriminus1*B3_aug-C_aug*Apoweriminus2*B3_aug;
            CApoweriminus1B31 = C_aug*Apoweriminus1*B3_aug;
 end
        
%===================For Hd===================
if ndi_aug~=0    
            CApoweriminus1Bd = C_aug*Apoweriminus1*Bd_aug-C_aug*Apoweriminus2*Bd_aug;
            CApoweriminus1Bd1 = C_aug*Apoweriminus1*Bd_aug;
end
    end
    
H1_dummy((i-1)*ny_aug+1:i*ny_aug,1:nu_aug)= CApoweriminus1B1;
%===================For H11===================
H11((i-1)*ny_aug+1:i*ny_aug,1:nu_aug)= CApoweriminus1B11;
    
if nd_aug~=0
        H2_dummy((i-1)*ny_aug+1:i*ny_aug,1:nd_aug)= CApoweriminus1B2;
%===================For H21===================        
        H21((i-1)*ny_aug+1:i*ny_aug,1:nd_aug)= CApoweriminus1B21;
end

if nz_aug~=0
        H3_dummy((i-1)*ny_aug+1:i*ny_aug,1:nz_aug)= CApoweriminus1B3;
%===================For H31===================
        H31((i-1)*ny_aug+1:i*ny_aug,1:nz_aug)= CApoweriminus1B31;
end

if ndi_aug~=0
        Hd_dummy((i-1)*ny_aug+1:i*ny_aug,1:ndi_aug)= CApoweriminus1Bd;
%===================For Hd1===================        
        Hd1((i-1)*ny_aug+1:i*ny_aug,1:ndi_aug)= CApoweriminus1Bd1;
    end
end


% %%%============================== Generation of H1 Matrix =============================================================%%%

for i=1:m
    %=========================H1===========================================
    H1((i-1)*ny_aug+1:p*ny_aug,(i-1)*nu_aug+1:i*nu_aug)= H1_dummy(1:(p+1-i)*ny_aug,1:nu_aug);
    
end
% ========================Change of last col of H1=========================

Apoweriminus1=eye(size(A_aug));
for i =1:p-m
    Apoweriminus1 = Apoweriminus1*A_aug;
    CApoweriminus1B1 = C_aug*Apoweriminus1*B1_aug;
    H1((m+i-1)*ny_aug+1:(m+i)*ny_aug,(m-1)*nu_aug+1:m*nu_aug)=CApoweriminus1B1;
    
end

% %%%============================== Generation of H2, H3, Hd Matrix =============================================================%%%

for i=1:p

    %=========================H2==================================
    H2((i-1)*ny_aug+1:p*ny_aug,(i-1)*nd_aug+1:i*nd_aug)=H2_dummy(1:(p+1-i)*ny_aug,1:nd_aug);
    
    %=========================H3==================================
    
    H3((i-1)*ny_aug+1:p*ny_aug,(i-1)*nz_aug+1:i*nz_aug)=H3_dummy(1:(p+1-i)*ny_aug,1:nz_aug);
    
    %=========================Hd==================================
    Hd((i-1)*ny_aug+1:p*ny_aug,(i-1)*ndi_aug+1:i*ndi_aug)=Hd_dummy(1:(p+1-i)*ny_aug,1:ndi_aug);
    
end 


%% ============================================== Generation of E1bar, E2bar, E3bar, E4bar, E5bar, Edbar in Equations (A.5-A.8)=========================================================
E1_bar= zeros(p*ne,m*nu_aug);
E2_bar= zeros(p*ne,p*nd_aug);
E3_bar = zeros(p*ne,p*nz_aug);
E4_bar = zeros(p*ne,p*ny_aug);
E5_bar = zeros(p*ne,1);
Ed_bar = zeros(p*ne,p*ndi_aug);
for i=1:p
    if i<=m
        E1_bar((i-1)*ne+1:i*ne,(i-1)*nu_aug+1:i*nu_aug)=-MldModel.E1;
    else
        E1_bar((i-1)*ne+1:i*ne,(m-1)*nu_aug+1:m*nu_aug)=-MldModel.E1;
    end
    E2_bar((i-1)*ne+1:i*ne,(i-1)*nd_aug+1:i*nd_aug)= MldModel.E2;
    E3_bar((i-1)*ne+1:i*ne,(i-1)*nz_aug+1:i*nz_aug)= MldModel.E3;
    E4_bar((i-1)*ne+1:i*ne,(i-1)*ny_aug+1:i*ny_aug)= -MldModel.E4;
    E5_bar((i-1)*ne+1:i*ne,1)= MldModel.E5;
    Ed_bar((i-1)*ne+1:i*ne,(i-1)*ndi_aug+1:i*ndi_aug)= MldModel.Ed;
end

if anticipation==0
Hd = Hd1;
Ed_bar = [MldModel.Ed; zeros((p*size(MldModel.Ed,1)-size(MldModel.Ed,1)),ndi_aug)]; 
Hd_bar=[zeros(ny_aug,ndi_aug);Hd(1:(p-1)*ny_aug,:)];
elseif anticipation==1
    Hd_bar=[zeros(ny_aug,p*ndi_aug);Hd(1:(p-1)*ny_aug,:)];
end




%% ===========================Generation of H1_bar, H2_bar, H3_bar, Hd_bar, H11_bar, H21_bar, H31_bar H_d1 bar and phi_bar in Equations (A.13-A.14)==========================================

H1_bar=[zeros(ny_aug,m*nu_aug);H1(1:(p-1)*ny_aug,:)];
H2_bar=[zeros(ny_aug,p*nd_aug);H2(1:(p-1)*ny_aug,:)];
H3_bar=[zeros(ny_aug,p*nz_aug);H3(1:(p-1)*ny_aug,:)];

H11_bar=[zeros(ny_aug,nu_aug);H11(1:(p-1)*ny_aug,:)];
H21_bar=[zeros(ny_aug,nd_aug);H21(1:(p-1)*ny_aug,:)];
H31_bar=[zeros(ny_aug,nz_aug);H31(1:(p-1)*ny_aug,:)];
Hd1_bar=[zeros(ny_aug,ndi_aug);Hd1(1:(p-1)*ny_aug,:)];

phi_bar=[C_aug; phi(1:(p-1)*ny_aug,:)];         % See Equation (A.14)


%% ======Generation of Epsilon_1,Epsilon_2, Epsilon_3, Epsilon_4, Epsilon_5, Epsilon_d, Epsilon_41, Epsilon_42, Epsilon_43, Epsilon_4d in Equations (A.9-A.12)==============================

Epsilon_1 =(E4_bar*H1_bar+E1_bar);
Epsilon_2 =(E4_bar*H2_bar+E2_bar);
Epsilon_3 =(E4_bar*H3_bar+E3_bar);
Epsilon_4=(E4_bar*phi_bar);
Epsilon_5 = E5_bar;
Epsilon_d = (E4_bar*Hd_bar+Ed_bar);

Epsilon_41 = E4_bar*H11_bar;
Epsilon_42 = E4_bar*H21_bar;
Epsilon_43 = E4_bar*H31_bar;
Epsilon_4d = E4_bar*Hd1_bar;

%=====================Ru and Ru0 matrices=========================
Rudummy=zeros(m*nu_aug,nu_aug);
ru=eye(nu_aug);
Rudummy(1:nu_aug,:)=ru;
Rudummy(nu_aug+1:2*nu_aug,:)=-ru;
Ru=zeros(m*nu_aug,m*nu_aug);
for i=1:m
    %=========================Ru===============================
    Ru((i-1)*nu_aug+1:m*nu_aug,(i-1)*nu_aug+1:i*nu_aug)=Rudummy(1:(m+1-i)*nu_aug,1:nu_aug);
end
%=========================Ru0==============================
Ru0=zeros(m*nu_aug,nu_aug);
Ru0(1:nu_aug,:)=ru;

a11=H1'*Wy*H1;
a12=H1'*Wy*H2;
a13=H1'*Wy*H3;
a21=H2'*Wy*H1;
a22=H2'*Wy*H2;
a23=H2'*Wy*H3;

a31=H3'*Wy*H1;
a32=H3'*Wy*H2;
a33=H3'*Wy*H3;

if slack_var==0
Qhat =[a11+Wu+(Ru'*Wdu*Ru), a12, a13;a21, a22+Wd, a23;a31, a32, a33+Wz]; % Wdu is weight on deltaU ans Wu is weight on U(general)
elseif slack_var==1    
 a14 = zeros(size(a13,1),size(Q_slack,2));
 a24 = zeros(size(a23,1),size(Q_slack,2));
 a34 = zeros(size(a33,1),size(Q_slack,2));
 a41 = a14';
 a42 = a24';
 a43 = a34';
 a44 = Q_slack; 
Qhat =[a11+Wu+(Ru'*Wdu*Ru), a12, a13, a14;a21, a22+Wd, a23, a24;a31, a32, a33+Wz, a34;a41, a42, a43, a44]; % Wdu is weight on deltaU ans Wu is weight on U(general)
end
% ========================--(General)=====================================

if slack_var==0
S1=[Epsilon_1 Epsilon_2 Epsilon_3];   %system constraint of MLD
S2=[H1 H2 H3];                        % corresponding to o/p constraint
S3=[eye(m*size(B1_aug,2)) zeros(m*size(B1_aug,2),p*size(B2_aug,2)) zeros(m*size(B1_aug,2),p*size(B3_aug,2))]; % corresponding to  constraint on I/P
S4=[Ru zeros(m*size(B1_aug,2),p*size(B2_aug,2)) zeros(m*size(B1_aug,2),p*size(B3_aug,2))];% corresponding to constraint on DELTA_U

Sb=[S1;S2;-S2;S3;-S3;S4;-S4];
elseif slack_var==1
  S1=[Epsilon_1, Epsilon_2, Epsilon_3, zeros(size(Epsilon_1,1),size(Q_slack,2))]; %system constraint of MLD 
  S2=[H1, H2, H3 -eye(size(H1,1),size(Q_slack,2))]; 
  S3 = [eye(m*size(B1_aug,2)), zeros(m*size(B1_aug,2),p*size(B2_aug,2)), zeros(m*size(B1_aug,2),p*size(B3_aug,2)) zeros(m*size(B1_aug,2),size(Q_slack,2))];
  S4 = [Ru, zeros(m*size(B1_aug,2),p*size(B2_aug,2)), zeros(m*size(B1_aug,2),p*size(B3_aug,2)) zeros(size(Ru,1),size(Q_slack,2))]; % corresponding to constraint on DELTA_U
  S5 = -[H1, H2, H3 eye(size(H1,1),size(Q_slack,2))];
  S6 = [zeros(ny_aug*p,size(Epsilon_1,2)),zeros(ny_aug*p,size(Epsilon_2,2)),zeros(ny_aug*p,size(Epsilon_3,2)), -eye(ny_aug*p,size(Q_slack,2))];  
  Sb=[S1; S2; S5; S3;-S3; S4; -S4; S6];
end

end
