function [D,V,Q] = find_markov_data(E,length,lanes,v0,p0,time,f)

n = numel(E);

D = cell(4,10);
V = cell(4,10);
Q = cell(4,10);

for i = 1:4
    for j = 1:10
        if j == 1
            d0 = 0.01;
        else
            d0 = 0.1*(j-1);
        end
        d = d0.*ones(n,1);
        [D{i,j},V{i,j},Q{i,j}] = markov_chain(d,E,length,lanes,v0(i)...
            ,p0(i),time,f);
        i
        j
    end
end
        