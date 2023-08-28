function Cruies_measurse_Vector
% Cumulative number of confirmed cases over Jan 20-Feb 19 for different 
%percentages of the protection measures

format long
%-------------Initial Conditions-----
 N=3711;
  n=2;
 N_p=2666;%Total population of passengers
 N_c=1045; %Total population of crew members

  %%%%%% Initial Conditions for the first part

 E00=0*ones(n,1);
 A00=0*ones(n,1);
I00=[0;1];
 R00=0*ones(n,1);
 S00=[N_c;N_p]-(E00+A00+I00+R00);
IC_firstpart_Vector=reshape([S00';E00';A00';I00';R00'],[],1);

%%%%%%%% Initial conditions for the second part

E0=[0.51;5];
A0=[2;4];
I0=[1;3];
R0=0*ones(n,1);
S0=[N_c;N_p]-(E0+A0+I0+R0);
IC_secondpart_Vector=reshape([S0';E0';A0';I0';R0'],[],1);
 
%-----------------------------------
 
t0=linspace(0,16,17); % Time interval for the first part of the voyage
 t1=linspace(17,31,15); % Time interval for the second part of the voyage

  op = odeset('RelTol',1e-5, 'AbsTol',1e-6);
 %-------- Mask Measure-----------
 etam=0.5; % Mask wear
 etan=0; % No vaccination
  
 s=[ 0.10 0.20 0.30 0.50 0.70 0.9]; %Percentage of Mask Measure

 for i=1:length(s)

  param_measure=s(i);

  param=[etan,etam,param_measure];

  %%%%%% Cumulative confirmed cases for the first part

  
[t,x_1]=ode45(@(t,x_1)Cruies_firstpart_measure_Vector(t,x_1,param),t0, IC_firstpart_Vector,op); %Jan 20-Feb 4
cumcase_first_Vector=cumsum(x_1(:,3:5:5*n)+x_1(:,4:5:5*n));
first_1=zeros(15,2);
first_1(1,1)=cumcase_first_Vector(17,1);
first_1(1,2)=cumcase_first_Vector(17,2);
first_1;

%%%%%% Cumulative confirmed cases for the second part

 IC_secondpartnew=in_condition_Vector_Measure( param, IC_firstpart_Vector); % Calcutes the initial condtions for the second part of the voyage for different measure
[t,x_2]=ode45(@(t,x_2)Cruies_secondpart_measure_Vector(t,x_2,param),t1,IC_secondpartnew,op);%Feb 5-Feb 20
cumcase_mask_Vector=cumsum(x_2(:,3:5:5*n)+x_2(:,4:5:5*n)+first_1);
cumcase_mask_Vector_final(i)=cumcase_mask_Vector(15,1)+cumcase_mask_Vector(15,2);

 end  
%------ Booster Measure-----
 etam=0; % No mask wear
 etan=0.74; % Booster measure
  
 s=[0.10 0.20 0.30 0.50 0.70 0.9]; %percentage of Booster Measure
 for i=1:length(s)

  param_measure=s(i);

  param=[etan,etam,param_measure];

     %%%%%% Cumulative confirmed cases for the first part

  
[t,x_1]=ode45(@(t,x_1)Cruies_firstpart_measure_Vector(t,x_1,param),t0, IC_firstpart_Vector,op); %Jan 20-Feb 4
cumcase_first_Vector=cumsum(x_1(:,3:5:5*n)+x_1(:,4:5:5*n));
first_1=zeros(15,2);
first_1(1,1)=cumcase_first_Vector(17,1);
first_1(1,2)=cumcase_first_Vector(17,2);
first_1;

%%%%%% Cumulative confirmed cases for the second part


IC_secondpartnew=in_condition_Vector_Measure( param, IC_firstpart_Vector);
[t,x_2]=ode45(@(t,x_2)Cruies_secondpart_measure_Vector(t,x_2,param),t1,IC_secondpartnew,op);%Feb 5-Feb 20
cumcase_booster_Vector=cumsum(x_2(:,3:5:5*n)+x_2(:,4:5:5*n)+first_1);
cumcase_booster_Vector_final(i)=cumcase_booster_Vector(15,1)+cumcase_booster_Vector(15,2);
 end
 %------- Vaccine Measure------
 etam=0; % No mask
 etan=0.94; % Vaccine measure
  
 s=[ 0.10 0.20 0.30 0.50 0.70 0.9]; %percentage of Vaccine Measure

 for i=1:length(s)

  param_measure=s(i);

  param=[etan,etam,param_measure];

     %%%%%% Cumulative confirmed cases for the first part

   
[t,x_1]=ode45(@(t,x_1)Cruies_firstpart_measure_Vector(t,x_1,param),t0, IC_firstpart_Vector,op); %Jan 20-Feb 4
cumcase_first_Vector=cumsum(x_1(:,3:5:5*n)+x_1(:,4:5:5*n));
first_1=zeros(15,2);
first_1(1,1)=cumcase_first_Vector(17,1);
first_1(1,2)=cumcase_first_Vector(17,2);
first_1;

%%%%%% Cumulative confirmed cases for the second part

IC_secondpartnew=in_condition_Vector_Measure( param, IC_firstpart_Vector);
[t,x_2]=ode45(@(t,x_2)Cruies_secondpart_measure_Vector(t,x_2,param),t1,IC_secondpartnew,op);%Feb 5-Feb 20
cumcase_vaccine_Vector=cumsum(x_2(:,3:5:5*n)+x_2(:,4:5:5*n)+first_1);
cumcase_vaccine_Vector_final(i)=cumcase_vaccine_Vector(15,1)+cumcase_vaccine_Vector(15,2);
 end
%------------ Vaccine+ Mask Measure------

 etam=0.55; % Mask wear
 etan=0.94;% Vaccine measure
  
 s=[ 0.10 0.20 0.30 0.50 0.70 0.90];  %percentage of Vaccine+ Mask Measure

 for i=1:length(s)

  param_measure=s(i);

  param=[etan,etam,param_measure];

 %%%%%% Cumulative confirmed cases for the first part

     
 
[t,x_1]=ode45(@(t,x_1)Cruies_firstpart_measure_Vector(t,x_1,param),t0, IC_firstpart_Vector,op); %Jan 20-Feb 4
cumcase_first_Vector=cumsum(x_1(:,3:5:5*n)+x_1(:,4:5:5*n));
first_1=zeros(15,2);
first_1(1,1)=cumcase_first_Vector(17,1);
first_1(1,2)=cumcase_first_Vector(17,2);
first_1;

%%%%%% Cumulative confirmed cases for the second part
 IC_secondpartnew=in_condition_Vector_Measure( param, IC_firstpart_Vector);
[t,x_2]=ode45(@(t,x_2)Cruies_secondpart_measure_Vector(t,x_2,param),t1,IC_secondpartnew,op);%Feb 5-Feb 20
cumcase_vaccine_mask_Vector=cumsum(x_2(:,3:5:5*n)+x_2(:,4:5:5*n)+first_1);
cumcase_vaccine_mask_Vector_final(i)=cumcase_vaccine_mask_Vector(15,1)+cumcase_vaccine_mask_Vector(15,2);
 end
 %--------- Booster+Mask Measure-------
etam=0.55; % Mask wear
 etan=0.74; % Booster measure
  
 s=[ 0.10 0.20 0.30 0.50 0.70 0.90]; %percentage of mesk measures

 for i=1:length(s)

  param_measure=s(i);

  param=[etan,etam,param_measure];

   %%%%%% Cumulative confirmed cases for the first part
     
   
[t,x_1]=ode45(@(t,x_1)Cruies_firstpart_measure_Vector(t,x_1,param),t0, IC_firstpart_Vector,op); %Jan 20-Feb 4
cumcase_first_Vector=cumsum(x_1(:,3:5:5*n)+x_1(:,4:5:5*n));
first_1=zeros(15,2);
first_1(1,1)=cumcase_first_Vector(17,1);
first_1(1,2)=cumcase_first_Vector(17,2);
first_1;

%%%%%% Cumulative confirmed cases for the second part
IC_secondpartnew=in_condition_Vector_Measure( param, IC_firstpart_Vector);
[t,x_2]=ode45(@(t,x_2)Cruies_secondpart_measure_Vector(t,x_2,param),t1,IC_secondpartnew,op);%Feb 5-Feb 20
cumcase_booster_mask_Vector=cumsum(x_2(:,3:5:5*n)+x_2(:,4:5:5*n)+first_1);
cumcase_booster_mask_Vector_final(i)=cumcase_booster_mask_Vector(15,1)+cumcase_booster_mask_Vector(15,2);
 end

Y=[cumcase_mask_Vector_final;cumcase_booster_Vector_final;cumcase_booster_mask_Vector_final;cumcase_vaccine_Vector_final;...
   cumcase_vaccine_mask_Vector_final]
 bar (Y)

  set(gca, 'XTickLabel',{'Mask', 'Booster' , 'Booster+Mask','Vaccine','Vaccine+Mask' })
legend('10%','20%','30%','50%','70%','90%')
ylabel('Number of  Cumulative Confirmed Cases over Jan 20-Feb 19')
title('Protection Measures for both Parts of the Voyage ')