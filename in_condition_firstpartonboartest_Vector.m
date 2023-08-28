function y=in_condition_firstpartonboartest_Vector( param, IC_firstpart)
% Calculate the initial conditions for the first part after RA testing
% measure on the second day of the voyage
format long

xi_s=0.73; % Sensitivity  of  RA test for symptomatic individuals
xi_a=0.57; % Sensitivity  of  RA test for asymptomatic individuals

n=2;
t0_1=linspace(0,2,3);
op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[t,x]=ode45(@(t,x)Cruies_firstpartonboardtest_Vector(t,x,param),t0_1,IC_firstpart,op);
S_onbaordtest=[x(3,1);x(3,11)];
E_onbaordtest=[x(3,2);x(3,12)];
A1_onbaordtest=[(1-xi_a)*x(3,3);(1-xi_a)*x(3,13)];
I1_onbaordtest=[(1-xi_s)*x(3,4);(1-xi_s)*x(3,14)];
A2_onbaordtest=[x(3,5);x(3,15)];
I2_onbaordtest=[x(3,6);x(3,16)];
R_onbaordtest=[x(3,7);x(3,17)];
A_iso=[xi_a*x(3,3);xi_a*x(3,13)];
I_iso=[xi_s*x(3,4);xi_s*x(3,14)];
R_iso=[0;0];

y1=[S_onbaordtest'; E_onbaordtest' ;A1_onbaordtest';I1_onbaordtest';A2_onbaordtest';I2_onbaordtest';R_onbaordtest';...
    A_iso';I_iso';R_iso'];
y=reshape(y1,[],1);
end