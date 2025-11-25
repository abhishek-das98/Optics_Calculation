function plot_polarization_ellipse(J, nPoints)
    % Plots the polarization ellipse corresponding to a Jones vector J.
    % nPoints controls the time sampling used to trace the ellipse.

    % If number of points is not provided, use a default value
    if nargin < 2 || isempty(nPoints)
        nPoints = 1000;
    end

    % Ensure J is a 2x1 Jones vector
    J = J(:);
    if numel(J) ~= 2
        error('J must be a 2x1 Jones vector [E_x; E_y].');
    end

    % Remove global phase so that the first nonzero component is real
    if abs(J(1)) >= 1e-14
        phi0 = angle(J(1));
    else
        phi0 = angle(J(2));
    end
    J = J * exp(-1i * phi0);

    % Normalize the Jones vector amplitude
    I = sqrt(abs(J(1))^2 + abs(J(2))^2);
    if I > 0
        J = J / I;
    end

    % Sample one optical cycle (dimensionless time ωt)
    t = linspace(0, 2*pi, nPoints);

    % Compute the real electric-field components (decreasing phase convention)
    E_t = real(J * exp(-1i * t));
    Ex = E_t(1, :);
    Ey = E_t(2, :);

    % Plot the polarization ellipse
    holdstate = ishold;
    plot(Ex, Ey, 'LineWidth', 2); hold on;
    axis equal; grid on;
    xlabel('E_x (arb.)');
    ylabel('E_y (arb.)');
    title('Polarization Ellipse');

    % Mark the field vector at t = 0
    plot(Ex(1), Ey(1), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 6);

    % Compute a small arrow showing the time evolution direction
    idx = round(numel(t)/8);     % pick a point not too close to the endpoints
    dx = Ex(idx+1) - Ex(idx);
    dy = Ey(idx+1) - Ey(idx);

    % Scale the arrow relative to ellipse size
    L = sqrt((max(Ex)-min(Ex))^2 + (max(Ey)-min(Ey))^2);
    arrowscale = 0.08 * L;

    % Normalize arrow direction and scale it
    vnorm = sqrt(dx^2 + dy^2);
    dx = dx / vnorm * arrowscale;
    dy = dy / vnorm * arrowscale;

    % Draw the arrow
    quiver(Ex(idx), Ey(idx), dx, dy, 0, ...
        'MaxHeadSize', 2.5, 'LineWidth', 1.8, 'Color', 'k');

    if ~holdstate
        hold off;
    end
end
