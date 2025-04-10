function s = change_bidsdir(s, bidsdir)
    
for f = fieldnames(s).'

    s_each = s.(f{1});
    if ischar(s_each) && contains(s_each, 'bids_dataset')
        s.(f{1}) = [bidsdir, extractAfter(s_each, 'bids_dataset')];
    end
    if isstruct(s_each)
        s.(f{1}) = change_bidsdir(s_each, bidsdir);
    end

end

end