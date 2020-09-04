function ODEsystem=orbit(v0,Jm,~)

G=6.67408*10^-11; % Gravitational constant in [m^3 kg^-1 s^-2].
M=1.98847*10^30;  % Solar mass in [kg].


r=v0(1);
aux=v0(2);

eq1=aux;              % r'
eq2=Jm*r^-3-G*M*r^-2; % aux'
ODEsystem=[eq1;eq2];
end