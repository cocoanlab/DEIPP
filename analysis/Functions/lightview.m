function lightview(az, el)

view(az, el);
delete(findobj(gca, 'Type', 'Light', '-depth', Inf));
lightangle(az, el);

end