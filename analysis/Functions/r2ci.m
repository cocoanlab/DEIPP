function ci = r2ci(r, n)

r = r(:);
n = n(:);

z = atanh(r);
alpha = 0.05;
zalpha = NaN(size(n));
zalpha(n>3) = (-erfinv(alpha - 1)) .* sqrt(2) ./ sqrt(n(n>3) - 3);
% zalpha(n>3) = norminv(1 - alpha/2) ./ sqrt(n(n>3) - 3);
ci = [tanh(z-zalpha) r tanh(z+zalpha)];

end