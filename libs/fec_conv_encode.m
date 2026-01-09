function [info_bits, coded_bits, pad_bits, sym_idx] = fec_conv_encode(num_data_syms, cfg)
    k = log2(cfg.M);
    target_coded = num_data_syms * k;
    info_len = floor(target_coded * cfg.fec.rate);
    info_bits = randi([0 1], info_len, 1);
    coded_bits = convenc(info_bits, cfg.fec.trellis);
    pad_bits = target_coded - numel(coded_bits);
    if pad_bits > 0
        coded_bits = [coded_bits; zeros(pad_bits, 1)];
    elseif pad_bits < 0
        coded_bits = coded_bits(1:target_coded);
        pad_bits = 0;
    end
    coded_mat = reshape(coded_bits, k, []).';
    sym_idx = bi2de(coded_mat, 'left-msb').';
end