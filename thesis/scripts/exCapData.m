function exCapData
    clearvars;
    close all;
    format shorte;
    
     filename = 'GRM31MR71H105KA88.txt';
     modelOrder  = 2;
     regression(filename,modelOrder);
%     [freq, cData] = getData(filename);
% 
%     figure;
%     [AX] = plotyy(freq,abs(cData),freq,rad2deg(phase(cData)),'loglog','semilogx');
% 
%     title('GRM31MR71H105KA88 Plot');
%     xlabel('f (Hz)');
%     ylabel(AX(1),'Magnitude (\Omega)');
%     ylabel(AX(2),'Phase (Degrees)');
%     set(AX(1),'FontSize', 15);  
%     set(AX(2),'FontSize', 15);
end
