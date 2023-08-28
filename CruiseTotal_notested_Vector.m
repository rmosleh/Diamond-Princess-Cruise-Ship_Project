function CruiseTotal_notested_Vector
% Cumulative number of confirmed cases over Jan 20-Feb 19 without RA testing measure onboard+
%impacts of the protection measures 
format long


 tau=0.96; % Sensitivity of  RT-PCR test 
 xi=0.73; % Sensitivity  of  RA test for symptomatic individuals
 

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
S0=[N_c;N_p]-(E0++A0+I0+R0);
IC_secondpart_Vector=reshape([S0';E0';A0';I0';R0'],[],1);
 
 %-----------Estimation-Non-Pretraveling tested--------------
 
%%%%%% Cumulative confirmed cases for the first part without pre-traveling testing


 t0=linspace(0,16,17); % Time interval for the first part of the voyage
 t1=linspace(17,31,15); % Time interval for the second part of the voyage


op = odeset('RelTol',1e-5, 'AbsTol',1e-6);

s= [0    0   0.74  0.74  0.94  0.94 ]; %etan
s1=[0  0.55  0   0.55   0    0.55 ]; % etam

  for i=1:length(s)

      etan=s(i);
      etam=s1(i);

       paramet=[etan, etam ] ;


 %%%%%% Cumulative confirmed cases for the first part

[t,x_1]=ode45(@(t,x_1)Cruies_firstpartnontested_Vector(t,x_1,paramet),t0,IC_firstpart_Vector,op); %Jan 20-Feb 4 
cumcase_first_Vector=cumsum(x_1(:,3:5:5*n)+x_1(:,4:5:5*n));
first_1=zeros(15,2);
first_1(1,1)=cumcase_first_Vector(17,1);
first_1(1,2)=cumcase_first_Vector(17,2);
first_1;
first_final=(cumcase_first_Vector(:,1)+cumcase_first_Vector(:,2));

 %%%%%% Cumulative confirmed cases for the second part
 

IC_secondpartnew_Vector=in_condition_Vector( paramet, IC_firstpart_Vector); % Calcutes the initial condtions for the second part of the voyage for different measure
[t,x_2]=ode45(@(t,x_2)Cruies_isolationcluster_Vector(t,x_2,paramet),t1,IC_secondpartnew_Vector,op);

cumcase_baseline_Vector=cumsum(x_2(:,3:5:5*n)+x_2(:,4:5:5*n)+first_1);
cumcase_baseline_Vector_final(i)=cumcase_baseline_Vector(15,1)+cumcase_baseline_Vector(15,2);

  end

%----------Pretraveling-PCR test-------

%%%%%% Cumulative confirmed cases for the first part with  RT-PCR pre-traveling testing

%%%%%%%% Initial conditions for RT-PCR pre-traveling testing measure

I00_PCR=(1-tau)*I00; % Calculates the initial infectees with RT-PCR pre-travelingg testing
IC_PCR=[S00'; E00'; A00'; I00_PCR'; R00'];
IC_PCR_Vector=reshape(IC_PCR,[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 s= [0    0   0.74  0.74  0.94  0.94 ]; %etan
  s1=[0  0.55  0   0.55   0    0.55 ]; % etam 

  for i=1:length(s)

      etan=s(i);
      etam=s1(i);

      paramet=[etan, etam ] ;


    %%%%%% Cumulative confirmed cases for the first part

[t,x_1_PCR]=ode45(@(t,x_1_PCR)Cruies_firstpartnontested_Vector(t,x_1_PCR,paramet),t0,IC_PCR_Vector,op); %Jan 20-Feb 4
cumcase_first_Vector_PCR=cumsum(x_1_PCR(:,3:5:5*n)+x_1_PCR(:,4:5:5*n));
first_1_PCR=zeros(15,2);
first_1_PCR(1,1)=cumcase_first_Vector_PCR(17,1);
first_1_PCR(1,2)=cumcase_first_Vector_PCR(17,2);
first_1_PCR;
first_final_PCR=(cumcase_first_Vector_PCR(:,1)+cumcase_first_Vector_PCR(:,2));


%%%%%% Cumulative confirmed cases for the second part
 
IC_secondpartnew_Vector=in_condition_Vector( paramet, IC_PCR_Vector);
[t,x_2_PCR]=ode45(@(t,x_2_PCR)Cruies_isolationcluster_Vector(t,x_2_PCR,paramet),t1,IC_secondpartnew_Vector,op);
cumcase_PCR_Vector=cumsum(x_2_PCR(:,3:5:5*n)+x_2_PCR(:,4:5:5*n)+first_1_PCR);
cumcase_PCR_Vector_final(i)=cumcase_PCR_Vector(15,1)+cumcase_PCR_Vector(15,2);

end        

% -------Pretraveling-RA test------

%%%%%% Cumulative confirmed cases for the first part with  RA pre-traveling testing
%%%%%%%% Initial conditions for RA pre-traveling testing measure

I00_RA=(1-xi)*I00;  % Calculates the initial infectees with RA pre-travelingg testing

IC_RA=[S00'; E00'; A00'; I00_RA'; R00'];
IC_RA_Vector=reshape(IC_RA,[],1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s= [0    0   0.74  0.74  0.94  0.94 ];%etan
  s1=[0  0.55  0   0.55   0    0.55 ];% etam

  for i=1:length(s)

      etan=s(i);
      etam=s1(i);

      paramet=[etan, etam ] ;

 %%%%%% Cumulative confirmed cases for the first part


[t,x_1_RA]=ode45(@(t,x_1_RA)Cruies_firstpartnontested_Vector(t,x_1_RA,paramet),t0,IC_RA_Vector,op); %Jan 20-Feb 4
cumcase_first_Vector_RA=cumsum(x_1_RA(:,3:5:5*n)+x_1_RA(:,4:5:5*n));
first_1_RA=zeros(15,2);
first_1_RA(1,1)=cumcase_first_Vector_RA(17,1);
first_1_RA(1,2)=cumcase_first_Vector_RA(17,2);
first_1_RA;
first_final_RA=(cumcase_first_Vector_RA(:,1)+cumcase_first_Vector_RA(:,2));


%%%%%% Cumulative confirmed cases for the second part

 IC_secondpartnew_Vector_RA=in_condition_Vector( paramet, IC_RA_Vector);
[t,x_2_RA]=ode45(@(t,x_2_RA)Cruies_isolationcluster_Vector(t,x_2_RA,paramet),t1,IC_secondpartnew_Vector_RA,op);
cumcase_RA_Vector=cumsum(x_2_RA(:,3:5:5*n)+x_2_RA(:,4:5:5*n)+first_1_RA);
cumcase_RA_Vector_final(i)=cumcase_RA_Vector(15,1)+cumcase_RA_Vector(15,2);

end        

Y=[cumcase_baseline_Vector_final ; cumcase_RA_Vector_final; cumcase_PCR_Vector_final]
bar(Y)

set(gca, 'XTickLabel',{'Non Pre-traveling Test ', ' Pre-Traveling RA Test ', 'Pre-Traveling PCR Test ' })
legend('Baseline','Mask','Booster','Boostere+Mask','Vaccine','Vaccine+Mask')
ylabel('Number of  Cumulative Confirmed Cases over Jan 20-Feb 19')
title('No Testing Measure on Board')

