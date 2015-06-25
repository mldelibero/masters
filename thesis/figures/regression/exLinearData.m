% The scripts in this file pertain to the modeling.tex file.
close all;

n = 10;
x = linspace(0,99,n);
b0 = rand * 10 -5 .* ones(1,n);
b1 = rand * 5 .* ones(1,n);
noise = rand(1,n) * 60 - 30;

figure;
plot(x,b0 + x .* b1 + noise,'O');
title('Example Plot');
xlabel('X');
ylabel('Y');

set(gca,'FontSize', 15);    