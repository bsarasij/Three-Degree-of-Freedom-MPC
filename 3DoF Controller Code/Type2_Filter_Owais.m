function [filt,unfilt] = Type2_Filter_Owais(n_w,alpha_d,Df,Df0_val)


[nt,nd]=size(Df);
% alpha_d=alpha_d(1:nd-1); % to adjust the different size of meass disturbance vector
na=length(alpha_d);
unfilt = Df(:,2:nd);
Df0=Df0_val*ones(n_w+1,1);%zeros(n_w+1,1);


Ae_f=zeros(n_w+1,n_w+1);
Ae_f(1,:)=[1 [1:n_w]*(-1/([1:n_w]*[1:n_w]'))]*alpha_d;
Ae_f(3:end,2:end-1)=eye(n_w-1);
Be_f=zeros(n_w+1,1);
Be_f(1,1)=(1-alpha_d)+(sum(1:n_w)/([1:n_w]*[1:n_w]'))*alpha_d;
Be_f(2,1)=1;
Ce_f=zeros(1,n_w+1);
Ce_f(1)=1;
    
Dfil_all=zeros(nt,1);
for jj=1:nt  
    Df1=Ae_f*Df0+Be_f*Df(jj,2:nd)';
    Df0=Df1;
    Dfil_all(jj,:)= Ce_f*Df1;
end
filt = Dfil_all; 

