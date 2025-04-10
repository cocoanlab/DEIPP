function s = calc_res_old(Y, Yfit, ntrain)

s = struct;

for j = 1:size(Y, 2)
    
    if j > ntrain
        ma = mean(cat(1, Y{:, 1:ntrain}));
    else
        ma = mean(cat(1, Y{:, setdiff(1:ntrain, j)}));
    end
    
    for i = 1:size(Y, 1)
        [a, b] = deal(Y{i,j}, Yfit{i,j});
        if all(isnan(a) | isnan(b))
            [a, b] = deal(NaN);
        else
            [a, b] = deal(a(~isnan(a) & ~isnan(b)), b(~isnan(a) & ~isnan(b)));
        end
        s.run.mse(i,j) = mean((a-b).^2);
        s.run.rsq(i,j) = 1 - mean((a-b).^2) ./ mean((a-ma).^2);
        [s.run.r_r(i,j), s.run.r_p(i,j)] = corr(a, b);
    end
    
    [a, b] = deal(cat(1, Y{:,j}), cat(1, Yfit{:,j}));
    if all(isnan(a) | isnan(b))
        [a, b] = deal(NaN);
    else
        [a, b] = deal(a(~isnan(a) & ~isnan(b)), b(~isnan(a) & ~isnan(b)));
    end
    s.ses.mse(j) = mean((a-b).^2);
    s.ses.rsq(j) = 1 - mean((a-b).^2) ./ mean((a-ma).^2);
    [s.ses.r_r(j), s.ses.r_p(j)] = corr(a, b);
    
end

ma = mean(cat(1, Y{:, 1:ntrain}));

[a, b] = deal(cat(1, Y{:,1:ntrain}), cat(1, Yfit{:,1:ntrain}));
if all(isnan(a) | isnan(b))
    [a, b] = deal(NaN);
else
    [a, b] = deal(a(~isnan(a) & ~isnan(b)), b(~isnan(a) & ~isnan(b)));
end
s.train.mse = mean((a-b).^2);
s.train.rsq = 1 - mean((a-b).^2) ./ mean((a-ma).^2);
[s.train.r_r, s.train.r_p] = corr(a, b);

if size(Y,2) > ntrain
    [a, b] = deal(cat(1, Y{:,ntrain+1:end}), cat(1, Yfit{:,ntrain+1:end}));
    if all(isnan(a) | isnan(b))
        [a, b] = deal(NaN);
    else
        [a, b] = deal(a(~isnan(a) & ~isnan(b)), b(~isnan(a) & ~isnan(b)));
    end
    s.test.mse = mean((a-b).^2);
    s.test.rsq = 1 - mean((a-b).^2) ./ mean((a-ma).^2);
    [s.test.r_r, s.test.r_p] = corr(a, b);
end

[a, b] = deal(cat(1, Y{:}), cat(1, Yfit{:}));
if all(isnan(a) | isnan(b))
    [a, b] = deal(NaN);
else
    [a, b] = deal(a(~isnan(a) & ~isnan(b)), b(~isnan(a) & ~isnan(b)));
end
s.all.mse = mean((a-b).^2);
s.all.rsq = 1 - mean((a-b).^2) ./ mean((a-ma).^2);
[s.all.r_r, s.all.r_p] = corr(a, b);

end