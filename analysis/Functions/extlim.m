function el = extlim(x, efac)

el = [min(x(:)) max(x(:))];

if nargin > 1
    d_el = diff(el);
    el = el + d_el * efac * [-1 1];
end

end