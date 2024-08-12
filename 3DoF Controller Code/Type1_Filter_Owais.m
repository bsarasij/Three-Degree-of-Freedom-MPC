function [filt] = Type1_Filter_Owais(alpha_d,Df)

[nt,nd]=size(Df);
alpha_d=alpha_d(1:nd-1); % to adjust the different size of meass disturbance vector
na=length(alpha_d);
Df0=zeros(na,1);
for jj=1:nt  
    Df1=diag(alpha_d)*Df0+(eye(na)-diag(alpha_d))*Df(jj,2:nd)';
    Df0=Df1;
    Dfil_all(jj,:)=Df1';  
end
filt = Dfil_all; 

 
