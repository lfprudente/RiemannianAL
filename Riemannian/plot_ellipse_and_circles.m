function plot_ellipse_and_circles(a, b, nballs, centers, radii, circle_colors,printPP,name)
    % a: semi-major axis length along x-axis
    % b: semi-minor axis length along y-axis
    % nballs: number of circles
    % centers: (nballs x 2) array with centers of the circles
    % radii: radii of the circles

    % Create figure
    figure;
    hold on;
    axis equal;
    
    % Plot the ellipse
    theta = linspace(0, 2*pi, 1000);
    x_ellipse = a * cos(theta);
    y_ellipse = b * sin(theta);
    plot(x_ellipse, y_ellipse, 'k', 'LineWidth', 2);

    % Plot the circles
    for k = 1:nballs
        x_circle = centers(k, 1) + radii * cos(theta);
        y_circle = centers(k, 2) + radii * sin(theta);
        %plot(x_circle, y_circle, 'k', 'LineWidth', 2);

        fill(x_circle, y_circle, circle_colors, 'EdgeColor', 'k', 'LineWidth', 2);

    end
    axis off
    
    % Set axis limits
    xlim([- a - radii, a + radii]);
    ylim([- b - radii, b + radii]);

    set(gcf, 'PaperPosition', [0 0 5 4]); 

    if ( printPP )
        print(name,'-depsc')
        saveas(gcf,strcat(name,'.png'))
    end
    
    hold off;
end
