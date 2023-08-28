function Cruies_cum_Vector
%Number of cumulative confirmed cases for each cluster separatly
format long

%--------------- Initial Conditions-------------
 N=3711;
  n=2;
 N_p=2666;%Total population of passengers
 N_c=1045; %Total population of crew members

%%%%%%%%%%%%% Initial condition for the first part

 E00=0*ones(n,1);
 A00=0*ones(n,1);
I00=[0;1];
 R00=0*ones(n,1);
 S00=[N_c;N_p]-(E00+A00+I00+R00);
 IC_firstpart_Vector=reshape([S00';E00';A00';I00';R00'],[],1);

 %%%%%%%%%%%% Initial condition for the second part


 E0=[0.51;5];
A0=[2;4];
 I0=[1;3];
 R0=0*ones(n,1);
S0=[N_c;N_p]-(E0+A0+I0+R0);
IC_secondpart_Vector=reshape([S0';E0';A0';I0';R0'],[],1);

 
%-----------Cumulative confirmed cases for both clusters---------------

 t0=linspace(0,16,17); %Time interval for the first part
 t1=linspace(17,31,15); %Time interval for the second part

  op = odeset('RelTol',1e-5, 'AbsTol',1e-6);

   s= [0    0   0.74  0.74  0.94  0.94 ]; %etan
  s1=[0  0.55  0   0.55   0    0.55 ];   % etam

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

 IC_secondpartnew_Vector=in_condition_Vector( paramet, IC_firstpart_Vector);
[t,x_2]=ode45(@(t,x_2)Cruies_isolationcluster_Vector(t,x_2,paramet),t1,IC_secondpartnew_Vector,op);

cumcase_baseline_Vector=cumsum(x_2(:,3:5:5*n)+x_2(:,4:5:5*n)+first_1);
cum_total_crew(i)=cumcase_baseline_Vector(15,1)
cum_total_pass(i)=cumcase_baseline_Vector(15,2)
cum_total(i)=cumcase_baseline_Vector(15,1)+cumcase_baseline_Vector(15,2)


  end
 %------------ Cumulative confirmed cases for each  cluster separately(protection measure only for the crew)

 s= [0    0   0.74  0.74  0.94  0.94 ]; %etan
  s1=[0  0.55  0   0.55   0    0.55 ];  %etam

  for i=1:length(s)
      etan=s(i);
      etam=s1(i);
paramet=[etan, etam ] ;

%%%%%% Cumulative confirmed cases for the first part

[t,x_crew]=ode45(@(t,x_crew)Cruies_firstpartnontested_Vector_crew(t,x_crew,paramet),t0,IC_firstpart_Vector,op); %Jan 20-Feb 4
cumcase_first_Vector_crew=cumsum(x_crew(:,3:5:5*n)+x_crew(:,4:5:5*n));
first_crew=zeros(15,2);
first_crew(1,1)=cumcase_first_Vector_crew(17,1);
first_crew(1,2)=cumcase_first_Vector_crew(17,2);
first_crew;

%%%%%% Cumulative confirmed cases for the second part


IC_secondpartnew_Vector_crew=in_condition_Vector_crew( paramet, IC_firstpart_Vector);
[t,x_crew_2]=ode45(@(t,x_crew_2)Cruies_isolationcluster_Vector_crew(t,x_crew_2,paramet),t1,IC_secondpartnew_Vector_crew,op);
cumcase_crew_2_Vector=cumsum(x_crew_2(:,3:5:5*n)+x_crew_2(:,4:5:5*n)+first_crew);
cumcase_crew_2_final(i)=cumcase_crew_2_Vector(15,1)+cumcase_crew_2_Vector(15,2);
cum_total_crew_cluster1(i)=cumcase_crew_2_Vector(15,1);
cum_total_pass_cluster1(i)=cumcase_crew_2_Vector(15,2);

  end
  %--------------Cumulative confirmed cases for each  cluster separately(protection measure only for the passengers)

  s= [0    0   0.74  0.74  0.94  0.94 ]; %etan
  s1=[0  0.55  0   0.55   0    0.55 ];  % etam

  for i=1:length(s)
      etan=s(i);
      etam=s1(i);
paramet=[etan, etam ] ;

%%%%%% Cumulative confirmed cases for the first part

[t,x_pass]=ode45(@(t,x_pass)Cruies_firstpartnontested_Vector_pass(t,x_pass,paramet),t0,IC_firstpart_Vector,op); %Jan 20-Feb 4
cumcase_first_Vector_pass=cumsum(x_pass(:,3:5:5*n)+x_pass(:,4:5:5*n));
first_pass=zeros(15,2);
first_pass(1,1)=cumcase_first_Vector_crew(17,1);
first_pass(1,2)=cumcase_first_Vector_crew(17,2);
first_pass;

%%%%%% Cumulative confirmed cases for the second part

IC_secondpartnew_Vector_pass=in_condition_Vector_pass( paramet, IC_firstpart_Vector);
[t,x_pass_2]=ode45(@(t,x_pass_2)Cruies_isolationcluster_Vector_pass(t,x_pass_2,paramet),t1,IC_secondpartnew_Vector_pass,op);
cumcase_pass_2_Vector=cumsum(x_pass_2(:,3:5:5*n)+x_pass_2(:,4:5:5*n)+first_pass);
cumcase_pass_2_final(i)=cumcase_pass_2_Vector(15,1)+cumcase_pass_2_Vector(15,2);
cum_total_crew_cluster2(i)=cumcase_pass_2_Vector(15,1);
cum_total_pass_cluster2(i)=cumcase_pass_2_Vector(15,2);
  end
%---------------------------------
Y_cluster_crew=[cum_total_crew_cluster1;cum_total_pass_cluster1]
Y_cluster_pass=[cum_total_crew_cluster2;cum_total_pass_cluster2]
Y=[cum_total_crew; cum_total_pass]

figure(1)
catStrArray = {'Crew ','Passenger '};
catStrArray ={'Baseline','Mask','Booster','Booster+Mask','Vaccine', 'Vaccine+Mask'}
catArray = categorical(catStrArray);       
catArray = reordercats(catArray,catStrArray);
bar(catArray,Y_cluster_pass)
 nModel = size(Y_cluster_pass,1);
nCat = size(Y_cluster_pass,2);
xPosAmpl = 0.3682626-0.3298725*exp(-0.407004*(nModel-1)); % position amplitude
 xPosInc = 2*xPosAmpl/(nModel-1);
 modelNames = [];
 for idxModel=1:nModel
    bar_xPos = 1:nCat;
    if nModel~=1
        bar_xPos = bar_xPos-xPosAmpl+(idxModel-1)*xPosInc;
    end
    text(bar_xPos,Y_cluster_pass(idxModel,:),num2str(Y_cluster_pass(idxModel,:)',...
        '%0.0f'),'vert','bottom','horiz','center'); 
   % modelNames{idxModel}=sprintf('model%d',idxModel);
end

 legend('Crew', 'Passengers')
 title('The Protection Measures  Are Applied only for the Passengers')
 ylabel('Number of  Cumulative Confirmed Cases over Jan 20-Feb 19')
figure(2)

catStrArray = {'Crew ','Passenger '};
catStrArray ={'Baseline','Mask','Booster','Booster+Mask','Vaccine', 'Vaccine+Mask'}
catArray = categorical(catStrArray);       
catArray = reordercats(catArray,catStrArray);
bar(catArray,Y_cluster_crew)
 nModel = size(Y_cluster_crew,1);
nCat = size(Y_cluster_crew,2);
xPosAmpl = 0.3682626-0.3298725*exp(-0.407004*(nModel-1)); % position amplitude
 xPosInc = 2*xPosAmpl/(nModel-1);
 modelNames = [];
 for idxModel=1:nModel
    bar_xPos = 1:nCat;
    if nModel~=1
        bar_xPos = bar_xPos-xPosAmpl+(idxModel-1)*xPosInc;
    end
    text(bar_xPos,Y_cluster_crew(idxModel,:),num2str(Y_cluster_crew(idxModel,:)',...
        '%0.0f'),'vert','bottom','horiz','center'); 
   % modelNames{idxModel}=sprintf('model%d',idxModel);
end

 legend('Crew', 'Passengers')
 title('The Protection Measures  Are Applied only for the Crew Memebers')
 ylabel('Number of  Cumulative Confirmed Cases')
 
figure(3)
catStrArray = {'Crew ','Passenger '};
catStrArray ={'Baseline','Mask','Booster','Booster+Mask','Vaccine', 'Vaccine+Mask'}
catArray = categorical(catStrArray);       
catArray = reordercats(catArray,catStrArray);
bar(catArray,Y)
 nModel = size(Y,1);
nCat = size(Y,2);
xPosAmpl = 0.3682626-0.3298725*exp(-0.407004*(nModel-1)); % position amplitude
 xPosInc = 2*xPosAmpl/(nModel-1);
 modelNames = [];
 for idxModel=1:nModel
    bar_xPos = 1:nCat;
    if nModel~=1
        bar_xPos = bar_xPos-xPosAmpl+(idxModel-1)*xPosInc;
    end
    text(bar_xPos,Y(idxModel,:),num2str(Y(idxModel,:)',...
        '%0.0f'),'vert','bottom','horiz','center'); 
   % modelNames{idxModel}=sprintf('model%d',idxModel);
end

 legend('Crew', 'Passengers')
 title('The Protection Measures  Are Applied  for Both Clusters')
 ylabel('Number of  Cumulative Confirmed Cases over Jan 20-Feb 19' )
 
 Y2=[cum_total]
 figure(4)
%catStrArray = {'Crew ','Passenger '};
catStrArray ={'Baseline','Mask','Booster','Booster+Mask','Vaccine', 'Vaccine+Mask'}
catArray = categorical(catStrArray);       
catArray = reordercats(catArray,catStrArray);
bar(catArray,Y2)
 nModel = size(Y2,1);
nCat = size(Y2,2);
xPosAmpl = 0.3682626-0.3298725*exp(-0.407004*(nModel-1)); % position amplitude
 xPosInc = 2*xPosAmpl/(nModel-1);
 modelNames = [];
 for idxModel=1:nModel
    bar_xPos = 1:nCat;
    if nModel~=1
        bar_xPos = bar_xPos-xPosAmpl+(idxModel-1)*xPosInc;
    end
    text(bar_xPos,Y2(idxModel,:),num2str(Y2(idxModel,:)',...
        '%0.0f'),'vert','bottom','horiz','center'); 
   % modelNames{idxModel}=sprintf('model%d',idxModel);
end

 
 title('The Protection Measures  Are Applied  for Both Clusters')
 ylabel('Number of  Cumulative Confirmed Cases over Jan 20-Feb 19')
 

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
dIdt(5:5:5*n)= [nu2 0;0 nu1]*flip(I(5:5:5*n))+r*(I(3:5:5*n)+I(4:5:5*n))-[nu1 0;0 nu2]*I(5:5:5*n);;

end
%-------------Feb5-Feb 20------------------
function dIdt=cruise31(t,I)
     dIdt=zeros(10*n,1);
    B1=I(3:10:10*n)+I(4:10:10*n)+I(5:10:10*n)+I(6:10:10*n);
    B2=I(8:10:10*n)+I(9:10:10*n);
l2=beta2*(B1./(1+c2*B1)+theta*(B2./(1+c3*B2)));% Force of Infection function
      L2=l2*(1-etan)*(1-etam);

dIdt(1:10:10*n)=[nu2 0;0 nu1]*flip(I(1:10:10*n))-w2*I(1:10:10*n)-L2.*I(1:10:10*n)-[nu1 0;0 nu2]*I(1:10:10*n);
dIdt(2:10:10*n)=[nu2 0;0 nu1]*flip(I(2:10:10*n))+L2.*I(1:10:10*n)-epsilon1*I(2:10:10*n)-[nu1 0;0 nu2]*I(2:10:10*n);
dIdt(3:10:10*n)=[nu2 0;0 nu1]*flip(I(3:10:10*n))+(1-p)*epsilon1*I(2:10:10*n)-epsilon2*I(3:10:10*n)-tau*alpha*I(3:10:10*n)-[nu1 0;0 nu2]*I(3:10:10*n);
dIdt(4:10:10*n)=[nu2 0;0 nu1]*flip(I(4:10:10*n))+p*epsilon1*I(2:10:10*n)-epsilon2*I(4:10:10*n)-tau*alpha*I(4:10:10*n)-[nu1 0;0 nu2]*I(4:10:10*n);
dIdt(5:10:10*n)=[nu2 0;0 nu1]*flip(I(5:10:10*n))+epsilon2*I(3:10:10*n)-r*I(5:10:10*n)-[nu1 0;0 nu2]*I(5:10:10*n);
dIdt(6:10:10*n)=[nu2 0;0 nu1]*flip(I(6:10:10*n))+epsilon2*I(4:10:10*n)-r*I(6:10:10*n)-[nu1 0;0 nu2]*I(6:10:10*n);
dIdt(7:10:10*n)= [nu2 0;0 nu1]*flip(I(7:10:10*n))+r*(I(3:10:10*n)+I(4:10:10*n))-[nu1 0;0 nu2]*I(7:10:10*n);
dIdt(8:10:10*n)= [nu2 0;0 nu1]*flip(I(8:10:10*n))+tau*alpha*I(3:10:10*n)-r*I(8:10:10*n)-[nu1 0;0 nu2]*I(8:10:10*n);
dIdt(9:10:10*n)= [nu2 0;0 nu1]*flip(I(9:10:10*n))+tau*alpha*I(4:10:10*n)-r*I(9:10:10*n)-[nu1 0;0 nu2]*I(9:10:10*n);
dIdt(10:10:10*n)= [nu2 0;0 nu1]*flip(I(10:10:10*n))+r*(I(8:10:10*n)+I(9:10:10*n))-[nu1 0;0 nu2]*I(10:10:10*n);




     end
end