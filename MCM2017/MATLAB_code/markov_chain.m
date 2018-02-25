function [D,V,Q] = markov_chain(d,E,length,lanes,v0,p0,time,f)

n = numel(E);

vol = length.*lanes;

D = zeros(n,time./1000);
V = zeros(n,time./1000);
Q = zeros(n,time./1000);

total = sum(vol.*d);
alpha = (circshift(E,1)./E) - 1;

for i = 2:time
    q = flow(d,v0,p0,lanes);
    outflow = q.*circshift((1-d),-1); % minus d
    inflow = circshift(q,1).*(1-d); % minus d
    excess = alpha.*outflow;
    if mod(i-1,1000) == 0
        Q(:,(i-1)./1000) = outflow;
    end
    
    d = d + (inflow-outflow-excess)./vol; 
    
    d = d.*(total/sum(d.*vol));
    
    if mod(i,1000) == 0
        D(:,i./1000) = d;
        v = velocity(d,v0,p0);
        V(:,i./1000) = v;
    end
end

q = flow(d,v0,p0,lanes); 

Q(:,time./1000) = q.*circshift((1-d),-1);


end
