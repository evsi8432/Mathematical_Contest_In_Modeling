
for i = 1:10
    d = 0.5.*rand(15,1);
    [D_5_dec0,V_5_dec0,Q_5_dec0] = markov_chain(d,E_520,length_520,lanes_520_dec,v0(1),p0(1),100000,1);
    hold on
    plot(mm_520,D_5_dec0(:,10))
    title('Steady-State Densities of 10 Simulations on SR-520 (Inital Density ~U(0,0.5))')
    ylabel('Density')
    xlabel('Mile Marker')
end