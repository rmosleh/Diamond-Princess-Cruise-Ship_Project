
function dIdt=Cruies_isolationcluster_seperate_Vector(t,I)
% Protection mesures are not considered for the second part of the voyage


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
r=0.07;


%------ System of the second part

 dIdt=zeros(5*n,1);
   B1=I(3:5:5*n)+I(4:5:5*n);
  l2=beta2*(B1./(1+c2*B1));
     
      


dIdt(1:5:5*n)=[nu2 0;0 nu1]*flip(I(1:5:5*n))-w2*I(1:5:5*n)-l2.*I(1:5:5*n)-[nu1 0;0 nu2]*I(1:5:5*n);
dIdt(2:5:5*n)=[nu2 0;0 nu1]*flip(I(2:5:5*n))+l2.*I(1:5:5*n)-epsilon0*I(2:5:5*n)-[nu1 0;0 nu2]*I(2:5:5*n);
dIdt(3:5:5*n)=[nu2 0;0 nu1]*flip(I(3:5:5*n))+(1-p)*epsilon0*I(2:5:5*n)-r*I(3:5:5*n)-[nu1 0;0 nu2]*I(3:5:5*n);
dIdt(4:5:5*n)=[nu2 0;0 nu1]*flip(I(4:5:5*n))+p*epsilon0*I(2:5:5*n)-r*I(4:5:5*n)-[nu1 0;0 nu2]*I(4:5:5*n);
dIdt(5:5:5*n)= [nu2 0;0 nu1]*flip(I(5:5:5*n))+r*(I(3:5:5*n)+I(4:5:5*n))-[nu1 0;0 nu2]*I(5:5:5*n);