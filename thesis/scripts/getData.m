function [w, cData, rData, iData] = getData(filename)
% This function parses data from an LTSpice file into a format that can be used for regression analysis
    data  = csvread(filename);
    w     = data(:,1) * 2 * pi;
    rData = data(:,2);
    iData = data(:,3);
    cData = rData + 1i * iData;
end
