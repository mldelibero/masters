function [w, cdata] = getData(filename)
    data      = csvread(filename);
    w         = data(:,1) * 2 * pi;
    real_data = data(:,2);
    imag_data = data(:,3);
    
    cdata = real_data + 1i * imag_data;
end