function v = velocity(d,v0,p0)

v = (0.8).*v0.*log((p0 + 0.31)./d);
v = v.*(v < 7) + 7.*(v > 7);
    
end