function y=in_conditionsecondpart_Vector( param, IC_firstpart)

% Calculates the initial conditions for the second part of the voyage 
%in on-board testing stream

t0_1=linspace(2,17,16);

param_1=[0 0 1]; % etan=0, etam=0
n=2;
N_p=2666;%Total population of passengers
 N_c=1045; %Total population of crew members
%--------- Initial Condition Second Part------
E0=[0.51;5];
A0=[2;4];
I0=[1;3];
 R0=0*ones(n,1);
S0=[N_c;N_p]-(E0+A0+I0+R0);

%------------------------------------


op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[t,x_1]=ode45(@(t,x_1)Cruies_firstpartonboardtest_Vector(t,x_1,param_1),t0_1,IC_firstpart,op);
[t,x]=ode45(@(t,x)Cruies_firstpartonboardtest_Vector(t,x,param),t0_1,IC_firstpart,op);

S_1=S0(1)/x_1(16,1);
S_2=S0(2)/x_1(16,11);
E_1=E0(1)/x_1(16,2);
E_2=E0(2)/x_1(16,12);
A_1=A0(1)/x_1(16,5);
A_2=A0(2)/x_1(16,15);
I_1=I0(1)/x_1(16,6);
I_2=I0(2)/x_1(16,16);

S_secondpart=[x(16,1)*S_1;x(16,11)*S_2]; 
E_secondpart=[x(16,2)*E_1;x(16,12)*E_2]; 
A_secondpart=[x(16,5)*A_1;x(16,15)*A_2]; 
I_secondpart=[x(16,6)*I_1;x(16,16)*I_2]; 
R_secondpart=[0;0];

y1=[S_secondpart'; E_secondpart' ;A_secondpart';I_secondpart';R_secondpart'];
y=reshape(y1,[],1);
end