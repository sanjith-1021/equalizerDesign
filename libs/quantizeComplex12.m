function y = quantizeComplex12(x)
%QUANTIZECOMPLEX12 Quantize complex values to signed 12-bit (Q11) int16.
    scale = 2^11 - 1;
    xr = max(min(real(x), 1), -1);
    xi = max(min(imag(x), 1), -1);
    re = int16(round(xr * scale));
    im = int16(round(xi * scale));
    y = complex(re, im);
end
