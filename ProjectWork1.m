%% ASTROPHYSICS and COSMOLOGY - Project Work 1

% Numerical integration of an elliptic orbit, by A. Justo and J. J. Ruiz.


%% General constants of the Solar System

G  = 6.67408*10^-11;       % Gravitational constant in [m^3 kg^-1 s^-2].
M  = 1.98847*10^30;        % Solar mass in [kg].
AU = 149597870700;         % Astronomical Unit in [m].
yr = 365.25*24*3600;       % Year in [s].


%% Data of the orbiting body: Halley's Comet

aph = 35.082*AU;           % Aphelion in [m].
a   = 17.834*AU;           % Major semiaxis in [m].
e   = 0.96714;             % Eccentricity (adimensional).
m   = 2.2*10^14;           % Mass in [kg].
Jm  = G*M*a*(1-e^2);       % Parameter equal to J^2/m^2 in SI units.
Jtheoretical = m*sqrt(Jm); % Angular momentum in SI units.


%% Estimation of the period by means of 3rd Kepler Law

% Period expected in [yr].
Ttheoretical=sqrt(4*pi^2*a^3/(G*M))/yr;


%% Parameters for the ODE solver

% First instant of integration.
t0 = 0;

% Final instant of integration, corresponding to 10 complete orbits.
tN = 10*Ttheoretical*yr;

% Time step.
k  = tN/10^5;

% Initial conditions. Impose that, at t=t0, the comet is located at the
% aphelium, where dr(t)/dt = 0.
initcond = [aph; 0];


%% Integration and results of the integration

% Obtain r(t) and dr(t)/dt.
rk4out = rk4(@orbit,initcond,Jm,k,tN,t0);
t      = rk4out(1,:);
r      = rk4out(2,:);
rdot   = rk4out(3,:);

% Numerically calculated aphelion.
aphelion=max(r);

% Numerically calculated perihelion.
perihelion=min(r);


%% Computing theta as a function of time

% Integrate dtheta(t)/dt = J/(m*r^2). Do the integration as the sum of
% trapezoid areas.

fvec=r.^-2;

% Preallocating for speed.
theta=zeros(1,length(t));

for i=2:length(t)
    if i==2
    theta(i) = sqrt(Jm)*k*(0.5*fvec(1) + 0.5*fvec(i));
    else
    theta(i) = sqrt(Jm)*k*(0.5*fvec(1) + sum(fvec(2:i-1)) + 0.5*fvec(i));
    end
end


%% Plot r(t) and theta(t) as a function of time

figure(1)
set(gcf,'Position',[175 350 1000 500])

subplot(2,1,1)
plot(t/yr,r/AU,'b',t/yr,aphelion*ones(1,length(t))/AU,'--g',t/yr,...
    perihelion*ones(1,length(t))/AU,'--r')
title('1P/Halley'); xlabel('t [years]'); ylabel('r [AU]'); grid
legend('r(t)','Aphelion','Perihelion')

subplot(2,1,2)
plot(t/yr,theta*180/pi); grid
title('1P/Halley'); xlabel('t [years]'); ylabel('\theta [º]')


%% Obtain the corresponding cartesian coordinates

x      = r.*cos(theta);
y      = r.*sin(theta);


%% Calculate the energy, which must be conserved

% Compute E(t).
E = 1/2*m*rdot.^2+1/2*Jm*m*r.^(-2)-G*M*m./r;

% Compute the % of lost E over all the integration.
lostE = 100*abs(1-E(end)/E(1));

% Since tN corresponds to 10 orbits...
lostEperOrbit = lostE/10;

% Compute the theoretical E, which must be conserved.
Etheoretical = -G*M*m/(2*a);


%% Calculate the angular momentum, which must be conserved

% Compute dtheta(t)/dt from the expression r(theta).
thetadot = a*(1-e^2)./(e*r.^2.*sqrt(1-(1/e*a*(1-e^2)./r-1/e).^2)).*rdot;

% Compute J(t).
J = m*r.^2.*thetadot;


%% Calculate the geometrical orbital parameters

% Numerically calculated major semiaxis in [AU].
a_num  = 0.5*(aphelion+perihelion)/AU

% Numerically calculated minor semiaxis in [AU].
b_num  = max(abs(y))/AU

% Numerically calculated eccentricity (three different ways).
e_num1 = sqrt(1-(b_num/a_num)^2)
e_num2 = sqrt(1-4*aphelion*perihelion/(aphelion+perihelion)^2)
e_num3 = (G^2*M^2*m^4/Jtheoretical^4+...
    2*m*Etheoretical/Jtheoretical^2)^(1/2)*Jtheoretical^2/(G*M*m^2)


%% Plot the energy and the angular momentum as a function of time

figure(2)
set(gcf,'Position',[175 350 1000 500])

subplot(2,1,1)
plot(t/yr,E,t/yr,Etheoretical*ones(1,length(t)),'--r')
title('1P/Halley')
xlabel('t [years]'); ylabel('E [J]'); grid
legend('Numerical E','Theoretical E')

subplot(2,1,2)
plot(t/yr,abs(J),t/yr,Jtheoretical*ones(1,length(t)),'--r')
title('1P/Halley')
xlabel('t [years]'); ylabel('|J| [kg m^2 s^-^1]'); grid
legend('Numerical |J|', 'Theoretical |J|')


%% Checking 1st Kepler Law: Elliptic trajectory with the Sun at a focus

% Check that the Sun is located at one of the focus.

% The first focus is located at the origin.
F1=[0,0];

% The second focus is located, from the aphelion, at a distance equal
% than the perihelion.
F2=[aphelion-perihelion,0];

% Now compute the sum of the distances of each point of the orbit from
% F1 and F2.
dF1=sqrt(x.^2+y.^2);
dF2=sqrt((x-F2(1)*ones(1,length(x))).^2+y.^2);

% Numerically calculated sum of the distance to both focuses, in [AU].
sum_numerical  =mean(dF1+dF2)/AU
sum_theoretical=2*a/AU


%% Checking 2nd Kepler Law: Equal areas are swept in equal times

% Compute the area swept in a given time.
x2=[0, x(1:2000)];
y2=[0, y(1:2000)];
area1=polyarea(x2/AU,y2/AU)

% Compute the area swept in the same time.
x3=[0, x(2001:2000*2)];
y3=[0, y(2001:2000*2)];
area2=polyarea(x3/AU,y3/AU)


%% Checking 3rd Kepler Law. Period

% Find the time it takes to come back to the aphelion. To do that, find
% the peak of r(t) when it arrives to the aphelion, as plotted in Fig (1).
[pks,locs] = findpeaks(r/AU, 'MinPeakHeight', aphelion/AU-1);

% The period is:
T=t(locs(1))/yr


%% Plot the resulting orbit and swept areas

figure(3)
set(gcf,'Position',[175 350 1000 500])

% Check that the points of the orbit fulfill expression of r(theta)
thetavec=linspace(0,2*pi,200);
r_theta=a*(1-e^2)./(1+e*cos(thetavec));
x_theta=-r_theta.*cos(thetavec);
y_theta=r_theta.*sin(thetavec);

plot(x/AU,y/AU,'k',x_theta/AU,y_theta/AU,'.b','markersize',10); hold on
plot(0,0,'.r','markersize',20);
title('1P/Halley');
xlabel('x [AU]'); ylabel('y [AU]'); axis equal; hold on; grid 


fill(x2/AU,y2/AU,'g'); hold on; grid
fill(x3/AU,y3/AU,'r')

legend('Orbit (1st way)','Orbit (2nd way)','Sun (focus)'); grid; xlim([-5 40])
text(20,1,['Area 1 = ' num2str(area1) ' AU^2'])
text(15,3,['Area 2 = ' num2str(area2) ' AU^2'])

text(15,-1.5,['a = ' num2str(a/AU) ' AU'])
text(15,-2.5,['T = ' num2str(T) ' years'])