clear
clc
close all
tic
global m q

angs2au=1.889726;      % angstrom to atomic unit (bohr radius)
au2fs=0.02418884;      % atomic unit to femtosecond
mpvsme=1822.888486424682;   % ratio of atomic mass and electron mass
hatoev=  27.211396; % taken from ksu webpage/ twice of ionization potential of H


% constants

m=[ 12 12 12 12  12 12 1.008 1.008 1.008 1.008 1.008  1.008]*mpvsme;
m = m';
q=[1 1 1 1 1 1 1 1 1 1 1 1];
% q=q'; % not necessary
ke = [];
p_ces_dis =[];

cart =        [   -1.02085400   -0.95725900   -0.00000200;
                   0.31872700   -1.36257100    0.00009400;
                   1.33948700   -0.40542800   -0.00008000;
                   1.02077400    0.95734600    0.00001400;
                  -0.31861300    1.36259500    0.00007300;
                  -1.33952100    0.40531900   -0.00007100;
                  -1.81625300   -1.70345300   -0.00006700;
                   0.56692500   -2.42460500    0.00008500;
                   2.38336400   -0.72128300   -0.00016300;
                   1.81635800    1.70333200    0.00001000;
                  -0.56706100    2.42456400    0.00002600;
                  -2.38332800    0.72143200   -0.00006200];
   



x = cart(:,1)*angs2au;         % Cartesian coordinates of the objects
y = cart(:,2)*angs2au ;        
z = cart(:,3)*angs2au;             

x_cm = sum(m.*x)/(sum(m));
y_cm = sum(m.*y)/(sum(m));
z_cm = sum(m.*z)/(sum(m));

x =  x - x_cm;
y =  y - y_cm;
z =  z - z_cm;
cart_cm = [x'; y'; z'];
cart_cm = reshape(cart_cm,1,[]);

% velocity components
vx = x.*0;
vy = y.*0;
vz = z.*0;
vel_cm = [vx'; vy'; vz'];
vel_cm = reshape(vel_cm,1,[]);

inits=[cart_cm vel_cm]; % initial conditions for each mass -- all in 1 column, 6 rows for ode solvers
%%
    
tspan=0:413413.789168889/100:413413.789168889*100;
options = odeset('reltol',1e-13,'abstol',1e-12); % tolerance for ode solver
[t,r] = ode45(@F,tspan,inits,options);
    
 
pos_vel_final = r(end,:);
pos = zeros(length(pos_vel_final)/6,3);
vel = pos;


%%
natoms = length(q); %can improve here

for j =1:natoms
    pos(j,1) = pos_vel_final(j*3-2);
    pos(j,2) = pos_vel_final(j*3-2+1);
    pos(j,3) = pos_vel_final(j*3-2+2);
    
    vel(j,1) = pos_vel_final(j*3-2+natoms*3);
    vel(j,2) = pos_vel_final(j*3-2+1+natoms*3);
    vel(j,3) = pos_vel_final(j*3-2+2+natoms*3);
        
end

format shortG
ke_mat = hatoev*0.5*m.*((vel.^2));
ke_f=[sum(ke_mat') sum(sum(ke_mat'))]

mom = m.*vel;
p_ces =reshape(mom', [],3*natoms);

dlmwrite('p_ces.csv', p_ces);
dlmwrite('ke_f.csv',ke);

%%
        r_ij=NaN(natoms,natoms);
        q_ij=r_ij; 

            
            for j=1:natoms
                for k=j+1:natoms
                     r_ij(j,k)= ( (cart(j,1)-cart(k,1))^2+(cart(j,2)-cart(k,2))^2+(cart(j,3)-cart(k,3))^2 )^0.5; %in A because taken from cart which is in A
                     q_ij(j,k)=q(j)*q(k);
                end
            end
            
            r_ij_inv = 1./r_ij;
            r_ij_inv=reshape( r_ij_inv,[],natoms*natoms);
            r_ij_inv=r_ij_inv(~isnan(r_ij_inv));
            
            q_ij =reshape( q_ij,[],natoms*natoms); 
            q_ij=q_ij(~isnan(q_ij));
            
            pe=14.4*(sum(q_ij.*r_ij_inv)) 


toc            
%%
function [deqn]=F(t,r)

global m q
natoms=length(m);

ii=[];
jj=[];
for i=1:natoms
    j=1:natoms;
    j(i)=[];
   for k=1:length(j)
        ii=[ii,i];
        jj=[jj,j(k)];
    end
end
index_pair=[ii' jj'];

ii=1:3*natoms;
ii=reshape(ii,[],natoms);

[mm,nn]=size(index_pair); %nn remain unsed

jj=NaN(3,1);
kk=jj;
for i=1:mm
   jj(:,i)=ii(:,index_pair(i,1));
    kk(:,i)=ii(:,index_pair(i,2));
end
jj=reshape(jj,[],1);
kk=reshape(kk,[],1);
jjj=reshape([index_pair(:,1) index_pair(:,1) index_pair(:,1)]',[],1); %for mass and charge
kkk=reshape([index_pair(:,2) index_pair(:,2) index_pair(:,2)]',[],1); %for mass and charge
i_ces=[jj kk jjj kkk];


vdef=zeros(3*natoms,1);

for i = 1:3*natoms
    vdef(i)=r(i+3*natoms); %generates first set of differential equations for initial velocity
end

    adef=[];
    adef_satom =[];

for j=1 :(natoms-1)*3:(natoms-1)*natoms*3  % can save invidual output ina matrix and then add them
    adef_satom_x = 0;
    adef_satom_y = 0;
    adef_satom_z = 0;
    for k =1:(natoms-1)
        ind=(k*3-2)+j-1;
     
        adef_satom_xk = ((q(i_ces(ind+0,3))*q(i_ces(ind+0,4))/m(i_ces(ind+0,3)))/(((r(i_ces(ind,1))-r(i_ces(ind,2))).^2 + (r(i_ces(ind+1,1))-r(i_ces(ind+1,2))).^2 + (r(i_ces(ind+2,1))-r(i_ces(ind+2,2))).^2).^(3/2)))*(r(i_ces(ind+0,1))-r(i_ces(ind+0,2)));
        adef_satom_yk = ((q(i_ces(ind+1,3))*q(i_ces(ind+1,4))/m(i_ces(ind+1,3)))/(((r(i_ces(ind,1))-r(i_ces(ind,2))).^2 + (r(i_ces(ind+1,1))-r(i_ces(ind+1,2))).^2 + (r(i_ces(ind+2,1))-r(i_ces(ind+2,2))).^2).^(3/2)))*(r(i_ces(ind+1,1))-r(i_ces(ind+1,2)));
        adef_satom_zk = ((q(i_ces(ind+2,3))*q(i_ces(ind+2,4))/m(i_ces(ind+2,3)))/(((r(i_ces(ind,1))-r(i_ces(ind,2))).^2 + (r(i_ces(ind+1,1))-r(i_ces(ind+1,2))).^2 + (r(i_ces(ind+2,1))-r(i_ces(ind+2,2))).^2).^(3/2)))*(r(i_ces(ind+2,1))-r(i_ces(ind+2,2)));
        
        adef_satom_x = adef_satom_x+adef_satom_xk;
        adef_satom_y = adef_satom_y+adef_satom_yk;
        adef_satom_z = adef_satom_z+adef_satom_zk;
        adef_satom =[adef_satom_x;adef_satom_y;adef_satom_z];
       
    end
    
    adef = [adef;adef_satom];
    
end

deqn=[vdef;adef]; 

end
