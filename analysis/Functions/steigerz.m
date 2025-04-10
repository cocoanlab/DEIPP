function [z, p] = steigerz(r_jk, r_jh, r_kh, n)

% Adapted from R package 'cocor: Comparing Correlations'
% 2023. 09. 27, Jae-Joong Lee

rm = (r_jk + r_jh) / 2;
c = (r_kh * (1 - 2*rm^2) - 1/2 * rm^2 * (1 - 2*rm^2 - r_kh^2)) / (1 - rm^2)^2;
z = ((atanh(r_jk) - atanh(r_jh)) * (n - 3).^0.5) / (2 - 2*c).^0.5;
p = 2 * (1 - normcdf(abs(z)));

end