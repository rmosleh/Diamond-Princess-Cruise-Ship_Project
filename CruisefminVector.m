function [p_estimate fval history] = CruisefminVector
% Estimates the parameters with the least square method by using fminsearch
% function

format long

history=[];
options = optimset('OutputFcn', @myoutput);

param1=[.003 .07 0.01 268000.214 10.111   0.18 .0002   ]; % Initial guess for the parameters

% nu1,nu2 ,w2, c1 ,c2, beta1, beta2

 time=linspace(0,31,32); %Time interval

[p_estimate,fval]=fminsearchbnd(@(param)fmin_estimate(param,time), param1, [0 0 0 10000 10  0  0  ], ...
    [1 3 1 3000000 10000   1 1  ],options) % Estimation process

%-------- Plotting process

data_cruise=load('DataCum.txt');
[x1,x2,x3]=fmin_estimate(p_estimate,time);
t0=linspace(0,16,17);
 t=linspace(17,31,15);
plot(t0, x2,'--','LineWidth',2);
hold on
plot(t, x3,t,data_cruise,'k .', 'markersize', 12,'LineWidth',2);
legend('Jan 20-Feb 4','Feb 5-Feb 19','Actual Data' ,'location', 'northwest' )
 xlabel('Days  (Jan 20-Feb 19)')
 ylabel('Number of the cumulative confirmed Cases')

 %----------- R_0 calculations----

r=0.07;
R_0_1=p_estimate(6)/r
R_0_2=p_estimate(7)/r

%--------- History of esitmation process
function stop = myoutput(p_estimate,optimvalues,state);
        stop = false;
        if isequal(state,'iter')
          history = [history; p_estimate];
        end
end

%------------ Confidence Intervals--
N_Sample=size(history,1);
%history(100:N_Sample,7);
% 
 N_1=size(history(70:N_Sample,1),1);
 N_2=size(history(110:N_Sample,2),1);
  N_3=size(history(60:N_Sample,3),1);
  N_4=size(history(80:N_Sample,4),1);
  N_5=size(history(70:N_Sample,5),1);
 N_7=size(history(20:N_Sample,7),1);
 
 %% SD
SD_nu1=std(history(70:N_Sample,1));
 SD_nu2=std(history(110:N_Sample,2));
 SD_w2=std(history(60:N_Sample,3));
  SD_c1=std(history(80:N_Sample,4));
   SD_c2=std(history(70:N_Sample,5));
 SD_beta1=std(history(:,6));
 SD_beta2=std(history(100:N_Sample,7));

%% Means

 mu_nu1=mean(history(70:N_Sample,1));
 mu_nu2=mean(history(110:N_Sample,2));
  mu_w2=mean(history(60:N_Sample,3));
 mu_c1=mean(history(80:N_Sample,4));
   mu_c2=mean(history(70:N_Sample,5));

 mu_beta1=mean(history(:,6));
  mu_beta2=mean(history(100:N_Sample,7));

%%  %Confidence Intervals

CI_nu1_left=mu_nu1-1.96*(SD_nu1/sqrt(N_1))
 CI_nu1_right=mu_nu1+1.96*(SD_nu1/sqrt(N_1))
 CI_nu2_left=mu_nu2-1.96*(SD_nu2/sqrt(N_2))  
 CI_nu2_right=mu_nu2+1.96*(SD_nu2/sqrt(N_2))
CI_w2_left=mu_w2-1.96*(SD_w2/sqrt(N_3))  
 CI_w2_right=mu_w2+1.96*(SD_w2/sqrt(N_3))

CI_c1_left=mu_c1-1.96*(SD_c1/sqrt(N_4))  
CI_c1_right=mu_c1+1.96*(SD_c1/sqrt(N_4))

CI_c2_left=mu_c2-1.96*(SD_c2/sqrt(N_5))  
CI_c2_right=mu_c2+1.96*(SD_c2/sqrt(N_5))

CI_beta1_left=mu_beta1-1.96*(SD_beta1/sqrt(N_Sample))  
CI_beta1_right=mu_beta1+1.96*(SD_beta1/sqrt(N_Sample))

CI_beta2_left=mu_beta2-1.96*(SD_beta2/sqrt(N_7))  
CI_beta2_right=mu_beta2+1.96*(SD_beta2/sqrt(N_7))







%-----------------------------------------
function [Y,z01,z11]=fmin_estimate(param,time)
nu1=param(1);
 nu2=param(2);
w2=param(3);
c1=param(4);
c2= param(5);
beta1= param(6);
beta2=   param(7);
p=0.821;
 epsilon0=0.17;
 epsilon1=0.598;
 epsilon2=0.5;
 r=0.07;
 tau=0.96; 
 n=2;
 etan=0;
 etam=0;
 %-------------Initial Conditions----
 N=3711;
 N_p=2666;%Total population of Passengers
 N_c=1045; %Total population odf
 E00=0*ones(n,1);
 A00=0*ones(n,1);
% I00=1*ones(n,1);
I00=[0;1];
 R00=0*ones(n,1);
 S00=[N_c;N_p]-(E00+A00+I00+R00);
 I01=[S00';E00';A00';I00';R00'];
I0r=reshape([S00';E00';A00';I00';R00'],[],1);

 E0=[0.51;5];
%A_01=0*ones(n,1);
A0=[2;4];
%I_01=0*ones(n,1);
 I0=[1;3];
 R0=0*ones(n,1);
 
 S0=[N_c;N_p]-(E0+A0+I0+R0);

 
 I05=[S0';E0';A0';I0';R0'];
I0r1=reshape([S0';E0';A0';I0';R0'],[],1);
%-----------------------------------
  t0=linspace(0,16,17); 
 t=linspace(17,31,15);

 op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[T,x]=ode45(@cruise3,t0,I0r,op);

z0=cumsum(x(:,3:5:5*n)+x(:,4:5:5*n));
z2_0=zeros(15,2);
z2_0(1,1)=z0(17,1);
z2_0(1,2)=z0(17,2);
z2_0;
z01=(z0(:,1)+z0(:,2));
[T1,y]=ode45(@cruise31,t,I0r1,op);
%y(:,3:10:20)+y(:,4:10:20)
  %z1=y(:,3:10:20)+y(:,4:10:20)+z2_0;
  z1=cumsum(y(:,3:5:5*n)+y(:,4:5:5*n)+z2_0);
z11=z1(:,1)+z1(:,2);
% z3=z2(15,1)+z2(15,2);
  %Y1=(z1(15,1)-82)^2;

%Y3=(z1(15,2)-537)^2;
%Y3=((z1(15,2)+z1(15,1))-621)^2;
Y=(z1(15,1)-83)^2+(z1(15,2)-538)^2;

% plot(t, z1,'LineWidth',2);
%-------- Jan 20-Feb 4----------

function dIdt=cruise3(t,I)

dIdt=zeros(5*n,1);
l1=beta1*(I(3:5:5*n)+I(4:5:5*n))./(1+c1*(I(3:5:5*n)+I(4:5:5*n)));% Force of Infection function

      L1=l1*(1-etan)*(1-etam);

%   
dIdt(1:5:5*n)=[nu2 0;0 nu1]*flip(I(1:5:5*n))-L1.*I(1:5:5*n)-[nu1 0;0 nu2]*I(1:5:5*n);
dIdt(2:5:5*n)=[nu2 0;0 nu1]*flip(I(2:5:5*n))+L1.*I(1:5:5*n)-epsilon0*I(2:5:5*n)-[nu1 0;0 nu2]*I(2:5:5*n);
dIdt(3:5:5*n)=[nu2 0;0 nu1]*flip(I(3:5:5*n))+(1-p)*epsilon0*I(2:5:5*n)-r*I(3:5:5*n)-[nu1 0;0 nu2]*I(3:5:5*n);
dIdt(4:5:5*n)=[nu2 0;0 nu1]*flip(I(4:5:5*n))+p*epsilon0*I(2:5:5*n)-r*I(4:5:5*n)-[nu1 0;0 nu2]*I(4:5:5*n);
dIdt(5:5:5*n)= [nu2 0;0 nu1]*flip(I(5:5:5*n))+r*(I(3:5:5*n)+I(4:5:5*n))-[nu1 0;0 nu2]*I(5:5:5*n);
end
%-------------Feb5-Feb 20------------------
function dIdt=cruise31(t,I)
   
   
     dIdt=zeros(5*n,1);
   B2=I(3:5:5*n)+I(4:5:5*n);
 
  %  B2=I(8:10:10*n)+I(9:10:10*n);
l2=beta2*(B2./(1+c2*B2));
      L2=l2*(1-etan)*(1-etam);




dIdt(1:5:5*n)=[nu2 0;0 nu1]*flip(I(1:5:5*n))-w2*I(1:5:5*n)-L2.*I(1:5:5*n)-[nu1 0;0 nu2]*I(1:5:5*n);
dIdt(2:5:5*n)=[nu2 0;0 nu1]*flip(I(2:5:5*n))+L2.*I(1:5:5*n)-epsilon0*I(2:5:5*n)-[nu1 0;0 nu2]*I(2:5:5*n);
dIdt(3:5:5*n)=[nu2 0;0 nu1]*flip(I(3:5:5*n))+(1-p)*epsilon0*I(2:5:5*n)-r*I(3:5:5*n)-[nu1 0;0 nu2]*I(3:5:5*n);
dIdt(4:5:5*n)=[nu2 0;0 nu1]*flip(I(4:5:5*n))+p*epsilon0*I(2:5:5*n)-r*I(4:5:5*n)-[nu1 0;0 nu2]*I(4:5:5*n);
dIdt(5:5:5*n)= [nu2 0;0 nu1]*flip(I(5:5:5*n))+r*(I(3:5:5*n)+I(4:5:5*n))-[nu1 0;0 nu2]*I(5:5:5*n);



     end
     
end
end