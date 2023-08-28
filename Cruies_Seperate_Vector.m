function Cruies_Seperate_Vector
%Impacts of the protection measures for only one part of the voyage
%with no testing measure on board

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

%------------------------------------------

 t0=linspace(0,16,17); % Time interval for the first part of the voyage
 t1=linspace(17,31,15); % Time interval for the second part of the voyage

  op = odeset('RelTol',1e-5, 'AbsTol',1e-6);

  % ---------- Protection Measrues only for the First Part-----
  
s= [0    0   0.74  0.74  0.94  0.94 ];%etan
  s1=[0  0.55  0   0.55   0    0.55 ];% etam

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

%%%%%% Cumulative confirmed cases for the second part

 IC_secondpartnew=in_condition_Vector( paramet, IC_firstpart_Vector);
[t,x_2]=ode45(@(t,x_2)Cruies_isolationcluster_seperate_Vector(t,x_2),t1,IC_secondpartnew,op); %Feb 5-Feb 19
cumcase_baseline_Vector=cumsum(x_2(:,3:5:5*n)+x_2(:,4:5:5*n)+first_1);
cumcase_Vector_final1(i)=cumcase_baseline_Vector(15,1)+cumcase_baseline_Vector(15,2)
      
  end
%---------- Protection Measrues only for the Second Part------

s= [0    0   0.74  0.74  0.94  0.94 ];%etan
  s1=[0  0.55  0   0.55   0    0.55 ];% etam

  for i=1:length(s)
      etan=s(i);
      etam=s1(i);
       paramet=[etan, etam ] ;
       

 %%%%%% Cumulative confirmed cases for the first part


[t,x_1]=ode45(@(t,x_1)Cruies_firstpartnontested_seperate_Vector(t,x_1),t0,IC_firstpart_Vector,op); %Jan 20-Feb 4
cumcase_first_Vector=cumsum(x_1(:,3:5:5*n)+x_1(:,4:5:5*n));
first_1=zeros(15,2);
first_1(1,1)=cumcase_first_Vector(17,1);
first_1(1,2)=cumcase_first_Vector(17,2);
first_1;

%%%%%% Cumulative confirmed cases for the second part


[t,x_2]=ode45(@(t,x_2)Cruies_isolationcluster_Vector(t,x_2,paramet),t1,IC_secondpart_Vector,op); %Feb 5-Feb 19
cumcase_Vector2=cumsum(x_2(:,3:5:5*n)+x_2(:,4:5:5*n)+first_1);
cumcase_Vector_final2(i)=cumcase_Vector2(15,1)+cumcase_Vector2(15,2)
       
 
  end
  
  Y=[cumcase_Vector_final1;cumcase_Vector_final2]
  bar(Y)
  set(gca, 'XTickLabel',{'Jan 20-Feb4, Protection Measures', ' Feb 5-Feb 19, Protection Measures' })
legend('Baseline','Mask','Booster','Boostere+Mask','Vaccine','Vaccine+Mask')
ylabel('Number of  Cumulative Confirmed Cases over Jan 20-Feb 19')
title('Protection Measures for Only One Part of the Voyage ')
