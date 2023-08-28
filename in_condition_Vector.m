function y=in_condition_Vector( param, IC)
% Calculate the inital conditions for the second part
t0_1=linspace(0,17,18);
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);

param_1=[0 0]; % etan=0, etam=0
n=2;
N_p=2666;%Total population of passengers
 N_c=1045; %Total population of crew members
%--------- Initial Condition Second Part------
E0=[0.51;5];
%A_01=0*ones(n,1);
A0=[2;4];
%I_01=0*ones(n,1);
 I0=[1;3];
 R0=0*ones(n,1);
S0=[N_c;N_p]-(E0+A0+I0+R0);

%------------------------------------

[t,x_1]=ode45(@(t,x_1)Cruies_firstpartnontested_Vector(t,x_1,param_1),t0_1,IC,op);
[t,x]=ode45(@(t,x)Cruies_firstpartnontested_Vector(t,x,param),t0_1,IC,op);

S_1=S0(1)/x_1(18,1);
S_2=S0(2)/x_1(18,6);
E_1=E0(1)/x_1(18,2);
E_2=E0(2)/x_1(18,7);
A_1=A0(1)/x_1(18,3);
A_2=A0(2)/x_1(18,8);
I_1=I0(1)/x_1(18,4);
I_2=I0(2)/x_1(18,9);



S_secondpart=[x(18,1)*S_1;x(18,6)*S_2]; 
E_secondpart=[x(18,2)*E_1;x(18,7)*E_2] ; 
A_secondpart=[x(18,3)*A_1;x(18,8)*A_2]; 
I_secondpart=[x(18,4)*I_1;x(18,9)*I_2];
R_secondpart=[0;0];


y=reshape([S_secondpart'; E_secondpart' ;A_secondpart';I_secondpart';R_secondpart'],[],1);
end