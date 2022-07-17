clear
clc
close all
tic
global m

angs2au=1.889726;      % angstrom to atomic unit (bohr radius)
au2fs=0.02418884;      % atomic unit to femtosecond
mpvsme=1822.888486424682;   % ratio of atomic mass and electron mass
hatoev=  27.211396; % taken from ksu webpage/ twice of ionization potential of H



% constants

m=[ 12 12 12 12  12 12 1.008 1.008 1.008 1.008 1.008  1.008]*mpvsme;
m = m';

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
   



X = cart(:,1)*angs2au;         % Cartesian coordinates of the objects
Y = cart(:,2)*angs2au ;        
Z = cart(:,3)*angs2au;             

X_cm = sum(m.*X)/(sum(m));
Y_cm = sum(m.*Y)/(sum(m));
Z_cm = sum(m.*Z)/(sum(m));

X =  X - X_cm;
Y =  Y - Y_cm;
Z =  Z - Z_cm;



%%
% velocity components
vx = X.*0;
vy = Y.*0;
vz = Z.*0;

inits=[X Y Z vx vy vz]; % initial conditions for each mass -- all in 1 column, 6 rows for ode solvers

    
time_span=0:413413.789168889/100:413413.789168889*100;
options = odeset('reltol',1e-13,'abstol',1e-12); % tolerance for ode solver
[t,r] = ode45(@nbody,time_span,inits,options);
    
 
pos_vel_final = r(end,:);
pos = zeros(length(pos_vel_final)/6,3);
vel = pos;

natoms = length(pos_vel_final)/6; %can improve here

for j =1:natoms
    pos(j,1) = pos_vel_final(j);
    pos(j,2) = pos_vel_final(j+natoms);
    pos(j,3) = pos_vel_final(j+natoms*2);
    
    vel(j,1) = pos_vel_final(j+natoms*3);
    vel(j,2) = pos_vel_final(j+natoms*4);
    vel(j,3) = pos_vel_final(j+natoms*5);
        
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

            
            for j=1:natoms
                for k=j+1:natoms
                     r_ij(j,k)= ( (cart(j,1)-cart(k,1))^2+(cart(j,2)-cart(k,2))^2+(cart(j,3)-cart(k,3))^2 )^0.5; %in A because taken from cart which is in A
                end
            end
            
            r_ij_inv = 1./r_ij;
            r_ij_inv=reshape( r_ij_inv,[],natoms*natoms);
            r_ij_inv=r_ij_inv(~isnan(r_ij_inv));
            pe=14.4*(sum(r_ij_inv))


toc

%%
function [dx] = nbody(t,x)

%REF: https://www.mathworks.com/matlabcentral/fileexchange/75202-n-body-simulation-with-ode45?s_tid=FX_rc1_behav
global  m

N = length(m); % The number of bodies
dx = zeros(6*N,1); % for each particle we need 6 equations: 3 for position, 3 for velocity; 

% computing distances r between bodies
% x=[X Y Z vx vy vz]
X=x(1:N);
Y=x(N+1:2*N);
Z=x(2*N+1:3*N);

DX=distance(X); % external function below
DY=distance(Y);
DZ=distance(Z);

r=(DX.^2+DY.^2+DZ.^2).^0.5; % matrix of distances between all particles

%% accelerations
   M=m*ones(1,N); % square matrix of masses, each row is m(i) repeated along columns (for next line)
   a=(1./M)./r.^2; % accelerations of each mass stuck here
   
   % get rid of diagonal entries (Inf - forces on bodies on itself)
   a(1:N+1:N*N) = 0;
    
   ax =a.* (DX./ r);
   ay =a.* (DY./ r);
   az =a.* (DZ./ r);
   
   ax(1:N+1:N*N) = 0;
   ay(1:N+1:N*N) = 0;
   az(1:N+1:N*N) = 0;
   
   % net accel on mass m(i) -- adding all a's along column number (direction 2)
   a_sx=sum(ax,2)';
   a_sy=sum(ay,2)';
   a_sz=sum(az,2)';

%% velocities

% Filling in the velocities from IC's
for i = 1:3*N
    dx(i) = x(3*N+i);    % dx is [vx1 vx2...vxN vy1 vy2...vyN vz1 vz2...vzN ax1 ax2 ...axN ay...etc az...etc]
end

for i=1:N
    dx(3*N+i)=a_sx(i)';
    dx(4*N+i)=a_sy(i)';
    dx(5*N+i)=a_sz(i)';
end
end



function [dist]=distance(X)
% this function calculates distance between two coordinates in a given vector X
L=length(X);
dist=zeros([1 L]);

j=1:L;
for k=1:L
  dist(j(k),j)=X(j(k))-X(j);
end
 
end