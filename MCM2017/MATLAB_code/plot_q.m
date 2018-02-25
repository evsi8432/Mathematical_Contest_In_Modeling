function plot_q(x,Q,d,f1,f2,density)

clf

figure(f1);

clf

for k = 1:4
    
subplot(2,2,k)

for i = 1:4
    plot(x,1.*mean(Q{i,d(k)}(:,50:100),2))
    hold on
end

legend('no SD Cars','10% SD Cars','50% SD Cars','90% SD Cars');
xlabel('Milemarker');
ylabel('Flow (Cars/Sec)');
ylim([0,1.5]);
title(sprintf('Flow of Cars Along SR-520 (Initial Density = %s)',density(k,:)))

end 

figure(f2)
clf

for k = 1:4

    subplot(4,1,k)
    plot(1000./60./60.*(1:100),Q{1,d(k)}(1,:));
    hold on
    for j = 2:numel(Q{1,d(k)}(:,1))
        plot(1000./60./60.*(1:100),Q{1,d(k)}(j,:))
    end
    title('Flow of Each Highway Segment Over Time (Initial Density = 0.01)');
    xlabel('Hours After Beginning of Sim');
    ylabel('Flow (Cars/Sec)')
end
