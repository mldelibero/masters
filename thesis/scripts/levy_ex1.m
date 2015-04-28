function levy_ex1()
clearvars;
close all;
format shorte;

filename = 'levy_ex1.csv';
modelOrder  = 2;
regression(filename,modelOrder);
end