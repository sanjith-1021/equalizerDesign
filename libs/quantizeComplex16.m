function y = quantizeComplex16(x)
%QUANTIZECOMPLEX16 Quantize complex values to signed 16-bit with Q15 scale.
    scale = 2^15 - 1;
    xr = max(min(real(x), 1), -1);
    xi = max(min(imag(x), 1), -1);
    xr = double(int16(round(xr * scale))) / scale;
    xi = double(int16(round(xi * scale))) / scale;
    y = complex(xr, xi);
end
