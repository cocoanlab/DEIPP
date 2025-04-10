function printprog(i, tot_i)

msg = sprintf(['%.' num2str(ceil(log10(tot_i))) 'd / %.d ...\n'], i, tot_i);
if i > 1; msg = [repmat('\b', 1, length(msg)) msg]; end
fprintf(msg);

end