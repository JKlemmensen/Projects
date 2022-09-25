clc
clear all
close all
%%
clc
L = 20;                     %Length of potential
global Interval
Interval = 0.01;             %Length of intervals

xarr = -L/2:Interval:L/2;   %Starting points
K = Kinetic(xarr);          %Kinetic energy

%% Harmonic oscillator
V = Potentialharmonic(xarr);%Potential energy
H = K+V;                    %Hamiltonian

[vec,val] = eig(H);         %Eigenfunctions, eigenvalues

%% Plotting the eigenfunctions
figure(1)
hold on
for i =1:2000             
plot(xarr,vec(:,i))
end
hold off

%% Finding the maximum value of $|\psi|^2$ as a function of n (eigenvalue index)
N=2100
[Y,I] = max(abs(vec(:,1:N)));       %I is index of the maximum value.
figure(3)
plot(1:N,xarr(I).^2,'o')

%% Probability of finding the particle in a position larger than X
X = 1;
N=50;
ProbabilityArray = abs(xarr) > X*ones(1,length(xarr));  %1 if x>X, otherwise 0.

Probability(1:N) = (ProbabilityArray*vec(:,1:N).^2./sum((vec(:,1:N).^2)))

plot(1:N,Probability)
title(['Probability of x>' num2str(X)])
ylabel('Probability')
xlabel('Index of eigenstate')



%% Gaussian Potential
clc
V = PotentialGaussian(10,xarr);%Potential energy
H = K+V;                    %Hamiltonian

[vec,val] = eig(H);         %Eigenfunctions, eigenvalues

%% Plotting the eigenfunctions
hold on
for i = 1:4
   plot(xarr,vec(:,i)); 
end
hold off

%% Minimum energy as a function of V0
V0 = 2:2:100
for i = V0
    VaryingPotential = PotentialGaussian(i,xarr);
    HVarying = K+VaryingPotential;
    [VecVarying,ValVarying] = eig(HVarying);
    MinimumEnergy(i/2) = ValVarying(1,1);        %Energy of ground state (E0)
    FirstEnergy(i/2) = ValVarying(2,2);

end
%%
GroundPotentialDifference = MinimumEnergy + V0;      %E0+V0
GroundFirstDifference = FirstEnergy - MinimumEnergy;      %E1-E0
%% Continued
figure(1)
plot(V0,MinimumEnergy,'o','markersize',8)
xlabel('Potentialet V_0')
ylabel('Grundtilstandsenergien E_0')
set(gca,'FontSize',20);
%%
figure(2)
hold on
plot(2:2:100,GroundPotentialDifference,'o','markersize',8)
%plot(2:2:100,floor(sqrt((2:2:100)/2)))
xlabel('Potentialet V_0')
ylabel('E_0+V_0')
%title('Plot of zero point energy E_0+V_0')
set(gca,'FontSize',20);
hold off

figure(3)
hold on
plot(2:2:100,GroundFirstDifference,'o','markersize',8)
%plot(2:2:100,floor(sqrt(2*(2:2:100))))
xlabel('Potentialet V_0')
ylabel('Energiforskellen E_1-E_0')
%title('Energy difference between ground and first excited state')
set(gca,'FontSize',20);
hold off
%% Scaling of bound states as V0 increases
j=1
for i = 1:2:400
    VaryingPotential = PotentialGaussian(i,xarr);
    HVarying = K+VaryingPotential;
    [VecVarying,ValVarying] = eig(HVarying);
    %MinimumEnergy(i/2) = ValVarying(1,1);        %Energy of ground state (E0)
    %FirstEnergy(i/2) = ValVarying(2,2);
    
    Scaling(j)= sum(sum(ValVarying<0));            %Number of bound states
    j=j+1;
end
%% Continued
%Scaling - sqrt(1:2:40)                   %Scaling contra sqrt(V0)
hold on
plot(1:2:400,Scaling,'markersize',8)
plot(1:2:400,sqrt(1.27*(1:2:400)),'markersize',8)
xlabel('Potentialet V_0')
%xlabel('$$\sqrt{f}$$','Interpreter','latex','FontSize',13)
ylabel('Antallet af bundne tilstande')
legend('Egentilstande','c\cdot sqrt(V_0)','Location','SouthEast')
hold off
set(gca,'FontSize',20);

%% Calculating photon wavelength when an electron jumps from 2->1

DeltaE = sqrt(1*(6.582e-16)^2/(2*511e3*(1e-9/(3e8))^2))

Lambda = 4.135e-16*3e8/DeltaE

%% Doing it numerically
clc
V = PotentialGaussian(10,xarr);%Potential energy
H = K+V;                    %Hamiltonian

[vec,val] = eig(H);         %Eigenfunctions, eigenvalues
%%
DeltaENum = val(2,2)-val(1,1)


%% Double well potential
clc
V = PotentialDouble(20,xarr);  %Potential energy
H = K+V;                    %Hamiltonian

[vec,val] = eig(H);         %Eigenfunctions, eigenvalues

%% Plotting the eigenfunctions
hold on
for i = 1:4
   plot(xarr,vec(:,i)); 
end
hold off

%% Potential exam

hbar = 6.582E-16 ; %In eV
m = 511E3;  %EV
sigma=0.4E-9;   %m
c=2.99E8

VeV = 4 %Measured in eV
E0 = m*sigma^2/(hbar*c)^2;
V0 = VeV*E0

V = Potential(V0,xarr);
%%
H = K+V;                    %Hamiltonian
[vec,val] = eig(H);         %Eigenfunctions, eigenvalues
%%
(val(2,2)-val(1,1))/E0

%% Calculating energy and wavelength for a hydrogen atom

hbar = 6.58E-16%*1.602E-19 ; %In eV
mElectron = 511E3%* 1.602E-19;  %eV
sigma=0.3E-9;   %m
c=3E8;
e = 1.602E-19;   %C
Epsilon = 8.854E-12*(1.602E-19);    %C^2/eV*m
n=3;
m=1E30;

EnergyHydrogen = mElectron*e^4/(2*(4*pi*Epsilon)^2*hbar^2*c^2)*(1/n^2-1/m^2)

lambda = hbar*2*pi*c/EnergyHydrogen

%% Infinite potential well

hbar = 6.58E-16; %In eV
mElectron = 511E3;  %EV
sigma=0.3E-9;   %m
c=3E8;
e = 1.602E-19;   %C
n = 1;
m = 0;
a = 0.8E-9;  %m

EnergyInfiniteWell = hbar^2*pi^2*c^2/(2*mElectron*a^2)* (n^2-m^2);

for n=2:3
    for m=1:n-1
        EnergyInfiniteWell = hbar^2*pi^2*c^2/(2*mElectron*a^2)* (n^2-m^2);
        lambda = hbar*2*pi*c/EnergyInfiniteWell
    end
end


%%

2*pi*hbar^3*51.1E9*(4*pi*Epsilon)^2*c^2/(mElectron*e^4)

