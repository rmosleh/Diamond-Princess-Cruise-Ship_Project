function CruiseTotal_onboardtested_Vector
% Cumulative number of confirmed cases over Jan 20-Feb 19 with RA testing measure on board+
%impacts of the protection measures


 tau=0.96; % Sensitivity of  RT-PCR test 
 xi_s=0.73; % Sensitivity  of  RA test for symptomatic individuals
 
 

 %-------------Initial Conditions-----
 N=3711;
  n=2;
 N_p=2666;%Total population of passengers
 N_c=1045; %Total population of crew members


 %%%%%% Initial Conditions for the first part 

 E00=0*ones(n,1);
 A00=0*ones(n,1);
I00=[0;0];
A00_1=0*ones(n,1);
I00_1=[0;1];
 R00=0*ones(n,1);
 A00_iso=0*ones(n,1);
 I00_iso=0*ones(n,1);
 R00_iso=0*ones(n,1);
 S00=[N_c;N_p]-(E00+A00_1+I00_1+A00+I00+R00+A00_iso+I00_iso+R00_iso);
IC_firstpart_Vector=reshape([S00';E00';A00';I00';A00_1';I00_1';R00';A00_iso';I00_iso';R00_iso'],[],1); % IC with RA testing on board 
IC_firstpart_Vector1=reshape([S00';E00';A00';I00_1';R00'],[],1); %IC without RA testing on board

%%%%%%%% Initial conditions for the second part

 E0=[0.51;5];
A0=[2;4];
 I0=[1;3];
 R0=0*ones(n,1);
S0=[N_c;N_p]-(E0+A0+I0+R0)
IC_secondpart_Vector=reshape([S0';E0';A0';I0';R0'],[],1);
 
 %-----------Baselines-Seperated--------------------

%%%%%% Cumulative confirmed cases for the first part without pre-traveling testing
%%%%%% With testing measure on board and without testing measure on board 

 t_00=linspace(0,16,17); % Time interval for the first part of the voyage
 t0=linspace(2,17,16); % Time interval for the first part of the voyage+ testing measure on the second day of the voyage
  t1=linspace(17,31,15); % Time interval for the second part of the voyage

  op = odeset('RelTol',1e-5, 'AbsTol',1e-6);

  s= [0    0   0.74  0.74  0.94  0.94 ]; %etan
  s1=[0  0.55  0   0.55   0    0.55 ]; % etam
 
  for i=1:length(s)

      etan=s(i);
      etam=s1(i);
      theta=1;

 paramet=[etan, etam, theta ] ;
 param=[etan, etam ] ;
 
    %%%%%% Cumulative confirmed cases for the first part with the testing
    %%%%%% measure on board


  IC_firstpartonboard_Vector=in_condition_firstpartonboartest_Vector( paramet, IC_firstpart_Vector); % IC for the first part after the isolation

  [t,x_1]=ode45(@(t,x_1)Cruies_firstpartonboardtest_Vector(t,x_1,paramet),t0,IC_firstpartonboard_Vector,op);%Jan 20-Feb 4
  cumcase_first_test=cumsum(x_1(:,3:10:10*n)+x_1(:,4:10:10*n)+x_1(:,5:10:10*n)+ ...
   x_1(:,6:10:10*n)+x_1(:,8:10:10*n)+x_1(:,9:10:10*n)) % With the isolation
  first_test=zeros(15,2);
 first_test(1,1)=cumcase_first_test(15,1);
 first_test(1,2)=cumcase_first_test(15,2);
 first_test;
  cumcase_first_test=cumcase_first_test(15,1)+cumcase_first_test(15,2);

   %%%%%% Cumulative confirmed cases for the second part with isolation
 
IC_secondpartnew_test=in_conditionsecondpart_Vector( paramet, IC_firstpartonboard_Vector); %IC for the second part with 
[t,x_2]=ode45(@(t,x_2)Cruies_isolationcluster_Vector(t,x_2,param),t1,IC_secondpartnew_test,op);
cumcase_baseline_Vector=cumsum(x_2(:,3:5:5*n)+x_2(:,4:5:5*n)+ first_test);
cumcase_baseline_Vector_final_test(i)=cumcase_baseline_Vector(15,1)+cumcase_baseline_Vector(15,2);

%%%%%% Cumulative confirmed cases for the first part without the testing measure on board
    


[t,x_3]=ode45(@(t,x_3)Cruies_firstpartnontested_Vector(t,x_3,param),t_00,IC_firstpart_Vector1,op); %Jan 20-Feb 4
cumcase_first_Vector=cumsum(x_3(:,3:5:5*n)+x_3(:,4:5:5*n)); %Without the isolation

first_1=zeros(15,2);
first_1(1,1)=cumcase_first_Vector(17,1);
first_1(1,2)=cumcase_first_Vector(17,2);
first_1;

%%%%%% Cumulative confirmed cases for the second part without isolation


IC_secondpartnew_Vector=in_condition_Vector( param, IC_firstpart_Vector1);
[t,x_2]=ode45(@(t,x_2)Cruies_isolationcluster_Vector(t,x_2,param),t1,IC_secondpartnew_Vector,op);

cumcase_baseline_Vector=cumsum(x_2(:,3:5:5*n)+x_2(:,4:5:5*n)+first_1);
cumcase_baseline_Vector_final_notest(i)=cumcase_baseline_Vector(15,1)+cumcase_baseline_Vector(15,2);
end
% 
%------------------ Ventilation effects-----------

s= [0    0   0.74  0.74  0.94  0.94 ]; %etan
  s1=[0  0.55  0   0.55   0    0.55 ]; % etam
  s2=[0 0.3 0.5 0.7 0.9 1]; %Ventilation effects
  for i=1:length(s)

      etan=0;
      etam=0;
      theta=s2(i);

 paramet=[etan, etam, theta ] ;
 param=[etan, etam ] ;
 
    %%%%%% Cumulative confirmed cases for the first part


  IC_firstpartonboard_Vector=in_condition_firstpartonboartest_Vector( paramet, IC_firstpart_Vector); % IC for the first part after the isolation

  [t,x_ven]=ode45(@(t,x_ven)Cruies_firstpartonboardtest_Vector(t,x_ven,paramet),t0,IC_firstpartonboard_Vector,op);%Jan 20-Feb 4
  cumcase_first_ven=cumsum(x_ven(:,3:10:10*n)+x_ven(:,4:10:10*n)+x_ven(:,5:10:10*n)+ ...
   x_ven(:,6:10:10*n)+x_ven(:,8:10:10*n)+x_ven(:,9:10:10*n)) % With the isolation


  first_ven=zeros(15,2);
first_ven(1,1)=cumcase_first_ven(15,1);
first_ven(1,2)=cumcase_first_ven(15,2);
first_ven;

   %%%%%% Cumulative confirmed cases for the second part
 
IC_secondpartnew_ven=in_conditionsecondpart_Vector( paramet, IC_firstpartonboard_Vector);
  [t,x_2_ven]=ode45(@(t,x_2_ven)Cruies_isolationcluster_Vector(t,x_2_ven,param),t1,IC_secondpartnew_ven,op);
cumcase_ven_Vector=cumsum(x_2_ven(:,3:5:5*n)+x_2_ven(:,4:5:5*n)+first_ven);
cumcase_ven_Vector_final_test(i)=cumcase_ven_Vector(15,1)+cumcase_ven_Vector(15,2);

  end


%  %-----------Estimation-Non-Pretraveling tested--------------




s= [0    0   0.74  0.74  0.94  0.94 ];%etan
  s1=[0  0.55  0   0.55   0    0.55 ];% etam

  for i=1:length(s)
      etan=s(i);
      etam=s1(i);
      theta=1;

 paramet=[etan, etam, theta ] ;
 param=[etan, etam ] ;

  %%%%%% Cumulative confirmed cases for the first part

 IC_firstpartonboard_Vector=in_condition_firstpartonboartest_Vector( paramet, IC_firstpart_Vector) ; 
  [t,x_test2]=ode45(@(t,x_test2)Cruies_firstpartonboardtest_Vector(t,x_test2,paramet),t0,IC_firstpartonboard_Vector,op);%Jan 20-Feb 4
  cumcase_first_test2=cumsum(x_test2(:,3:10:10*n)+x_test2(:,4:10:10*n)+x_test2(:,5:10:10*n)+ ...
      x_test2(:,6:10:10*n)+x_test2(:,8:10:10*n)+x_test2(:,9:10:10*n));

 first_test2=zeros(15,2);
first_test2(1,1)=cumcase_first_test2(15,1);
first_test2(1,2)=cumcase_first_test2(15,2);
first_test2;

 %%%%%% Cumulative confirmed cases for the second part
 
  IC_secondpartnew_test2=in_conditionsecondpart_Vector( paramet, IC_firstpartonboard_Vector);
  [t,x_test2]=ode45(@(t,x_test2)Cruies_isolationcluster_Vector(t,x_test2,paramet),t1,IC_secondpartnew_test2,op);
cumcase_baseline_Vector=cumsum(x_test2(:,3:5:5*n)+x_test2(:,4:5:5*n)+first_test2);
cumcase_baseline_Vector_final(i)=cumcase_baseline_Vector(15,1)+cumcase_baseline_Vector(15,2);

end        

% % %----------Pretraveling-PCR test-------

I00_PCR=(1-tau)*I00_1;  %IC with PCR pre-traveling test

IC_PCR=[S00'; E00'; A00'; I00'; A00_1'; I00_PCR'; R0'; A00_iso'; I00_iso'; R00_iso'];
IC_PCR_Vector=reshape(IC_PCR,[],1);

s= [0    0   0.74  0.74  0.94  0.94 ];%etan
  s1=[0  0.55  0   0.55   0    0.55 ];% etam

  for i=1:length(s)
      etan=s(i);
      etam=s1(i);
       theta=1;

 paramet=[etan, etam, theta ] ;
 param=[etan, etam ] ;


  %%%%%% Cumulative confirmed cases for the first part

 IC_firstpartonboard_Vector_PCR=in_condition_firstpartonboartest_Vector( paramet, IC_PCR_Vector) ; 
  [t,x_PCR]=ode45(@(t,x_PCR)Cruies_firstpartonboardtest_Vector(t,x_PCR,paramet),t0,IC_firstpartonboard_Vector_PCR,op);%Jan 20-Feb 4
  cumcase_first_PCR=cumsum(x_PCR(:,3:10:10*n)+x_PCR(:,4:10:10*n)+x_PCR(:,5:10:10*n)+ ...
      x_PCR(:,6:10:10*n)+x_PCR(:,8:10:10*n)+x_PCR(:,9:10:10*n));

 first_PCR=zeros(15,2);
first_PCR(1,1)=cumcase_first_PCR(15,1);
first_PCR(1,2)=cumcase_first_PCR(15,2);
first_PCR;
  IC_secondpartnew_PCR=in_conditionsecondpart_Vector( paramet, IC_firstpartonboard_Vector_PCR);
  [t,x_2_PCR]=ode45(@(t,x_2_PCR)Cruies_isolationcluster_Vector(t,x_2_PCR,param),t1,IC_secondpartnew_PCR,op);
cumcase_PCR_Vector=cumsum(x_2_PCR(:,3:5:5*n)+x_2_PCR(:,4:5:5*n)+first_PCR);
cumcase_PCR_Vector_final(i)=cumcase_PCR_Vector(15,1)+cumcase_PCR_Vector(15,2)

end        

% % %-------Pretraveling-RA test------

I00_RA=(1-xi_s)*I00_1;  %IC with PCR pre-traveling test


IC_RA=[S00'; E00'; A00'; I00'; A00_1'; I00_RA'; R0'; A00_iso'; I00_iso'; R00_iso'];
IC_RA_Vector=reshape(IC_RA,[],1);

s= [0    0   0.74  0.74  0.94  0.94 ];%etan
  s1=[0  0.55  0   0.55   0    0.55 ];% etam

  for i=1:length(s)
      etan=s(i);
      etam=s1(i);
       theta=1;

 paramet=[etan, etam, theta ] ;
 param=[etan, etam ] ;
 

 %%%%%% Cumulative confirmed cases for the first part

  IC_firstpartonboard_Vector_RA=in_condition_firstpartonboartest_Vector( paramet, IC_RA_Vector);      
[t,x_RA]=ode45(@(t,x_RA)Cruies_firstpartonboardtest_Vector(t,x_RA,paramet),t0,IC_RA_Vector,op);%Jan 20-Feb 4
  cumcase_first_RA=cumsum(x_RA(:,3:10:10*n)+x_RA(:,4:10:10*n)+x_RA(:,5:10:10*n)+ ...
      x_RA(:,6:10:10*n)+x_RA(:,8:10:10*n)+x_RA(:,9:10:10*n));


 first_RA=zeros(15,2);
first_RA(1,1)=cumcase_first_RA(15,1);
first_RA(1,2)=cumcase_first_RA(15,2);
first_RA;

%%%%%% Cumulative confirmed cases for the first part

  IC_secondpartnew_RA=in_conditionsecondpart_Vector( paramet, IC_firstpartonboard_Vector_RA);
  [t,x_2_RA]=ode45(@(t,x_2_RA)Cruies_isolationcluster_Vector(t,x_2_RA,param),t1,IC_secondpartnew_RA,op);
cumcase_RA_Vector=cumsum(x_2_RA(:,3:5:5*n)+x_2_RA(:,4:5:5*n)+first_RA);
cumcase_RA_Vector_final(i)=cumcase_RA_Vector(15,1)+cumcase_RA_Vector(15,2);
  end        

  %--------------Plotting 

Y_1=[cumcase_baseline_Vector_final_notest;cumcase_baseline_Vector_final_test]
figure(1)
bar(Y_1)
set(gca, 'XTickLabel',{'No Testing Measure on Board ','Testing Measure on Board' })
legend('Baseline','Mask','Booster','Boostere+Mask','Vaccine','Vaccine+Mask')
title(' Impacts of the Protection Measures without and with RA Testing Measure on Board')

%%%%%%%%%%%%%%%%%%%%%%%

Y_2=[cumcase_baseline_Vector_final;cumcase_RA_Vector_final;cumcase_PCR_Vector_final]
figure(2)
bar(Y_2)
set(gca, 'XTickLabel',{'Non Pre-Traveling Test ', ' Pre-Traveling RA Test ', 'Pre-Traveling PCR Test ' })
legend('Baseline','Mask','Booster','Boostere+Mask','Vaccine','Vaccine+Mask')
title(' Impacts of  the Protection Measures along with RA Testing Measure on Board ')


%%%%%%%%%%%%%%%%%

Y_3=[cumcase_ven_Vector_final_test]

figure(3)
bar(Y_3)
catStrArray ={'Baseline','\theta=0.3','\theta=0.5','\theta=0.7','\theta=0.9','\theta=1'}
catArray = categorical(catStrArray);       
catArray = reordercats(catArray,catStrArray);
bar(catArray,Y_3)
 nModel = size(Y_3,1);
nCat = size(Y_3,2);
xPosAmpl = 0.3682626-0.3298725*exp(-0.407004*(nModel-1)); % position amplitude
 xPosInc = 2*xPosAmpl/(nModel-1);
 modelNames = [];
 for idxModel=1:nModel
    bar_xPos = 1:nCat;
    if nModel~=1
        bar_xPos = bar_xPos-xPosAmpl+(idxModel-1)*xPosInc;
    end
    text(bar_xPos,Y_3(idxModel,:),num2str(Y_3(idxModel,:)',...
        '%0.0f'),'vert','bottom','horiz','center'); 
   % modelNames{idxModel}=sprintf('model%d',idxModel);
end

set(gca, 'XTickLabel',{'Baseline','\theta=0.3','\theta=0.5','\theta=0.7','\theta=0.9','\theta=1' })

ylabel('Number of  Cumulative Confirmed Cases over Jan 22-Feb 19')
title(' Impacts of Ventilation')
