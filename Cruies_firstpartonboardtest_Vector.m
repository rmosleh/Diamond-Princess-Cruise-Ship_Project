function dIdt=Cruies_firstpartonboardtest_Vector(t,I, param)
% This function solves the system for the first part of the voayge with 
% the  RA testing measure on board

%-------- Parameters for different measure
etan=param(1);
etam=param(2);
theta=param(3);
%----------- Parameters
n=2;
nu1=0.08104765372 ;
nu2=0.00028269383;
w2=0.03441628070;
c1=  10961.90480219778;
c2=  36.27769179422 ;
beta1= 0.22704791096;
beta2=  0.13753459741;
 p=0.821;
epsilon0=0.17;
 epsilon1=0.598;
 epsilon2=0.5;
 epsilon3=1;
 r=0.07;
 


N_p=2666;%Total population of passengers
 N_c=1045; %Total population of crew members
N=[N_c;N_p];

%--------------- System for the first part of the voyage
dIdt=zeros(10*n,1);

B1=I(3:10:10*n)+I(4:10:10*n)+I(5:10:10*n)+I(6:10:10*n);
       B2=I(8:10:10*n)+I(9:10:10*n);
       
        lambda_1=beta1*(B1./(1+c1*B1)+theta.*(B2./N));
         L1=lambda_1*(1-etan)*(1-etam); % Force of Infection function



dIdt(1:10:10*n)=[nu2 0;0 nu1]*flip(I(1:10:10*n))-L1.*I(1:10:10*n)-[nu1 0;0 nu2]*I(1:10:10*n);
dIdt(2:10:10*n)=[nu2 0;0 nu1]*flip(I(2:10:10*n))+L1.*I(1:10:10*n)-epsilon2*I(2:10:10*n)-[nu1 0;0 nu2]*I(2:10:10*n);
dIdt(3:10:10*n)=[nu2 0;0 nu1]*flip(I(3:10:10*n))+(1-p)*epsilon2*I(2:10:10*n)-epsilon3*I(3:10:10*n)-[nu1 0;0 nu2]*I(3:10:10*n);
dIdt(4:10:10*n)=[nu2 0;0 nu1]*flip(I(4:10:10*n))+p*epsilon2*I(2:10:10*n)-epsilon3*I(4:10:10*n)-[nu1 0;0 nu2]*I(4:10:10*n);
dIdt(5:10:10*n)= [nu2 0;0 nu1]*flip(I(5:10:10*n))+epsilon3*I(3:10:10*n)-r*I(5:10:10*n)-[nu1 0;0 nu2]*I(5:10:10*n);
dIdt(6:10:10*n)= [nu2 0;0 nu1]*flip(I(6:10:10*n))+epsilon3*I(4:10:10*n)-r*I(6:10:10*n)-[nu1 0;0 nu2]*I(6:10:10*n);
dIdt(7:10:10*n)= [nu2 0;0 nu1]*flip(I(7:10:10*n))+r*(I(5:10:10*n)+I(6:10:10*n))-[nu1 0;0 nu2]*I(7:10:10*n);
dIdt(8:10:10*n)= [nu2 0;0 nu1]*flip(I(8:10:10*n))-r*I(8:10:10*n)-[nu1 0;0 nu2]*I(8:10:10*n);
dIdt(9:10:10*n)= [nu2 0;0 nu1]*flip(I(8:10:10*n))-r*I(9:10:10*n)-[nu1 0;0 nu2]*I(9:10:10*n);
dIdt(10:10:10*n)= [nu2 0;0 nu1]*flip(I(8:10:10*n))+r*(I(8:10:10*n)+I(9:10:10*n))-[nu1 0;0 nu2]*I(10:10:10*n);