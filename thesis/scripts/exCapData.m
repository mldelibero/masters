function exCapData
    clearvars;
    close all;
    format shorte;

    filename = 'GRM31MR71H105KA88.txt';
    modelOrder  = 5;
    regression(filename,modelOrder);
end
