
%I skrev at det var tilladt at udvide som man lystede, så det har jeg gjort
%med et meget realistisk og livstruende eksempel

%Hvis opgaven skal løses ordentligt skal massen af overraskelsen blot
%sættes lig nul samt delete-funktionen slettes, således at planetbanernes
%spor ikke slettes aht. overskuelighed.

clear all
close all
clc
clear variables

Np = 7;         %Antal himmellegemer

% [Solen, Jorden, Mars, Merkur, stor planet 1, mørkt stof, stor planet 2];

m = [1, 3.004 * 10^(-6), 0.107 * 3.004 * 10^(-6), 0.815 * 3.004 * 10^(-6)...
    0.00001, 5, 0.0001];    %Alle masserne er opgivet i solmasser
%Hvis opgaven skal løses efter bogen skal 5 ændres til 0.


G = 4*(3.1415)^2;           %Gravitationskonstanten i AU, solmasser og år

vx = [0, 6.328, 0 , 0, 4, 50, -3];      %Hastighederne for planeterne
vy = [0, 0, 5.063 , 7.383,0, 5, 0, -6];


x = [0, 0, 1.524, 0.7233, 0, -70, 3];   %Deres startposition
y = [0, 1.000, 0, 0, 3, -4, 4];

l=1;                                    %Blot en tællevariabel

dt = 0.01;                              %Lille tidsskridt

figure(1)
xlabel('x [AU]','fontsize',14)
ylabel('y [AU]','fontsize',14)
title('Planetsimulering','fontsize',18)
hold on
plot(0,0,'or')
axis([-8 4 -4 8])
for j=0:dt:10
for i=1:Np
    ax = 0;
    ay = 0;
    
    for k=1:Np        %Evaluerer gravitationskraften fra alle objekterne
        
    if i ~= k         %Dog ikke sin egen
    r(k) = power((x(i)-x(k))^2 + (y(i)-y(k))^2,3/2);
    rx(k) = (x(k)-x(i));
    ry(k) = (y(k)-y(i));    %Nyttige variable
    
    
    %Finder summen af de accelerationer, som hver planet får fra de
    %omliggende
    ax = ax + (G * m(k))/(r(k)) * rx(k);
    ay = ay + (G * m(k))/(r(k)) * ry(k);
    end
    end
    
    %Sædvanlig udregning af afstand og hastighed i tidsrummet dt
    vx(i) = vx(i) + ax*dt;
    vy(i) = vy(i) + ay*dt;
    x(i) = x(i) + vx(i)*dt;
    y(i) = y(i) + vy(i)*dt;
    
    %Plotter punktet som findes ved hvert loop
    h(i,l) = scatter(x(i),y(i),5);
   
    %Sletter punktet der ligger 50 skridt tilbage, så plottet bliver mere
    %overskueligt. Kan slettes hvis planetbanernes form vil beholdes.
    if l > 50
        delete(h(i,l-50));
    end

end
pause(0.001)    %Pause
l = l+1;          %Tællevariablen stiger med 1


%Tekst til plottet
if j == 0.1
a = -7.8;
b = -3.6;
str = 'Wait for it...';
q1 = text(a,b,str,'FontSize',12);
end

if j == 0.8
a = -6;
b = -3.6;
str = 'Wait for it...';
q2 = text(a,b,str,'FontSize',12);
end

if j == 1.2
a = -4;
b = -3.6;
str = 'En kæmpe klump mørkt stof på fem solmasser!';
text(a,b,str,'FontSize',12);
end

end
hold off

