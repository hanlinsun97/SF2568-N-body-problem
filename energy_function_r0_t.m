clear
close all

% Position matrix
R_TT =[];
% Time array
TT = [];

% Constants
eps0 = 1;
rho0 = 1;
m = 1;
beta = (1/3)*rho0/(eps0*m);

% Integration step
h = 0.01;

% initial sphere
dr = 0.1;
r0 = dr:dr:1;

% Looping on time
for T=0:0.1:5

% R_T is the position at time T for different r0
R_T = zeros(1,length(r0));
% For every initial positions
for i=1:length(r0)
    t = 0;
    r = r0(i)+h;
    while(t<T)
        integrand = sqrt((r0(i)*r)/(2*beta*r0(i)^3*(r-r0(i))));
        t = t + integrand*h;
        r = r+h;
    end
    R_T(i) = r;
end
R_TT = [R_TT ; R_T];
TT = [TT t];
end

for i=1:length(r0)
    
    plot(TT,R_TT(:,i));
    hold on
    r0(i)

end
xlabel('time')
ylabel('position r')
pause


% SO NOW WE HAVE POSITION WITH T : we need to compute energy. 

% Velocity
V_TT = zeros(size(R_TT,1),size(R_TT,2));
for i=1:length(r0)
    V_TT(:,i) = sqrt(2*beta*r0(i)^3*(1/r0(i) - 1./R_TT(:,i)));
end

E_TT = 0.5 * m * V_TT.^2;

% Energy plot
figure

for i=1:length(r0)
    plot(TT,E_TT(:,i));
    
    hold on
end
xlabel('time')
ylabel('Energy')
pause


% FINALLY we want energy function of number of particles at time t1 (index)
t1 = 50;

N = size(r0);
figure

for i=1:length(r0)
    N(i) = 4*pi*(r0(i)^2 * dr - r0(i)*dr^2)*rho0 ;
    plot(E_TT(t1,i),N(i),'o')
    hold on
end
ylabel('Number of particle')
xlabel('Energy')
