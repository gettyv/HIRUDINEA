function sys_fig = graph_sys(N, n_v)
    % Create graph
    sys_fig = figure;
    set(sys_fig, 'Visible', 'off');
    hold on;
    scatter(N(1,:), N(2,:));
    
    % Plot members of N
    for i = 1:3:size(N, 2)-2
        line([N(1,i) N(1,i+1)], [N(2,i) N(2,i+1)], 'Color', 'b');
        line([N(1,i+1) N(1,i+2)], [N(2,i+1) N(2,i+2)], 'Color', 'b');
        line([N(1,i+2) N(1,i)], [N(2,i+2) N(2,i)], 'Color', 'b');
    end
    
    % Plot strings of N
    for i = 1:3:3*(n_v-1)
        line([N(1,i+4) N(1,i+1)], [N(2,i+4) N(2,i+1)], 'Color', 'r');
        line([N(1,i+3) N(1,i)], [N(2,i+3) N(2,i)], 'Color', 'r');
        line([N(1,i+4) N(1,i+2)], [N(2,i+4) N(2,i+2)], 'Color', 'r');
        line([N(1,i+3) N(1,i+2)], [N(2,i+3) N(2,i+2)], 'Color', 'r');
    end
    hold off;
    
    
    axis([-1 5 -3 3]);
    axis square;
    xlabel('X Position');
    ylabel('Y Position');
    title(append('System defined by N, where n_v = ',num2str(n_v)));
end

