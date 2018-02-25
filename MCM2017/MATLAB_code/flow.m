function q = flow(d,v0,p0,l)

q = (0.8).*l.*v0.*d.*log((p0 + 0.31)./(d));

end