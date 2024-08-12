clear
close all
clc
data=zeros(25,2);
data(2,:)=[8.5 7.5];
data(15,:)=[19 20];
data(18,:)=[1 50];
data(25,:)=[5 10];
Tsamp=60;
for iter=1:2
    [Flowrates,Controller_name]=HMPC_main_bioreactor(data(:,iter),Tsamp);
end

