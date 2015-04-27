function [freq, cdata] = getData(filename)
    data      = csvread(filename);
    freq      = data(:,1);
    real_data = data(:,2);
    imag_data = data(:,3);
    
    cdata = real_data + 1i * imag_data;
end