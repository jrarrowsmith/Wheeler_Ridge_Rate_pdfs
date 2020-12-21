function plot_connections(older_age,younger_age,older_pos,younger_pos, portion_of_connections, n, lineweight)
%Plot subset of connections between samples of two terraces

% for i=1:portion_of_connections:n %draw a portion of the connections
%     dt = older_age(i)-younger_age(i)
%     dE = older_pos(i)-younger_pos(i)
%     rate=dE/dt
%     if rate<=0 %we don't want negative rates
%         plot([younger_age(i) older_age(i)],[younger_pos(i) older_pos(i)], 'k-', 'LineWidth',lineweight)
%         drawnow
%     end
% end

for i=1:portion_of_connections:n %draw a portion of the connections
    dt = older_age(i)-younger_age(i);
    dE = older_pos(i)-younger_pos(i);
    rate=dE/dt;
	if rate>=0 %we don't want negative rates
        plot([younger_age(i) older_age(i)],[younger_pos(i) older_pos(i)], 'k-', 'LineWidth',lineweight)
        drawnow
    end
end

