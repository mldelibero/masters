% The scripts in this file pertain to the modeling.tex file.
close all;
clear all;

R = 10; 
L = 0.02; 
C = 1*10^-6;
n = 100;
f = logspace(1,6,n);
w = 2*pi*f;

Z = R + 1i * w * L + 1 ./ (1i * w * C);

figure;
[AX,H1,H2] = plotyy(f,abs(Z),f,rad2deg(phase(Z)),'loglog','semilogx');

title('RLC Plot');
xlabel('f (Hz)');
ylabel(AX(1),'Mag (\Omega)');
ylabel(AX(2),'Phase (degrees)');
set(AX(1),'FontSize', 15);  
set(AX(2),'FontSize', 15);  