function res = ynu_trace_chief_arbitrary(surfaces, u_hat, stop_idx, n_air, z_obj)
% Traces the chief ray in Y–ν (height–optical-direction) form
% for an arbitrary stop position in a rotationally symmetric system.
%
% surfaces: [R, t_to_next, n_after]  (N×3)
% u_hat   : chief ray angle in object space (rad) in air
% stop_idx: index of the surface that acts as the stop (1 <= stop_idx <= N)
% n_air  : refractive index in object/object-space (typically 1)
% z_obj  : axial distance from object plane to surface 1 (in medium n_air)

% Chief ray at the stop:
y_stop  = 0;                 % by definition: passes through stop center
nu_stop = n_air * u_hat;     % ν = n * u; stop medium is assumed air here

N       = size(surfaces, 1);
R       = surfaces(:, 1);
t       = surfaces(:, 2);
n_after = surfaces(:, 3);

% Pre-compute n_before and power φ for each surface
n_before = zeros(N, 1);
phi      = zeros(N, 1);
n_curr   = n_air;
for i = 1:N
    n_before(i) = n_curr;
    if isinf(R(i))
        phi(i) = 0;  % plane or stop: no power
    else
        phi(i) = (n_after(i) - n_curr) / R(i);
    end
    n_curr = n_after(i);
end

% Storage for ray data
y_before_mm = nan(N, 1);
nu_before   = nan(N, 1);
u_before_rad = nan(N, 1);
nu_after    = nan(N, 1);
u_after_rad = nan(N, 1);

% ---------------------------
% Backward trace from stop to earlier surfaces
% ---------------------------
y  = y_stop;
nu = nu_stop;

for i = (stop_idx-1):-1:1
    % Inverse translation through medium n_after(i)
    y  = y - (nu / n_after(i)) * t(i);
    % Inverse refraction at surface i
    nu = nu + y * phi(i);

    % State just before surface i (in forward sense)
    y_before_mm(i) = y;
    nu_before(i)   = nu;
    u_before_rad(i)= nu / n_before(i);

    % State just after surface i (for completeness)
    nu_after(i)    = nu - y * phi(i);
    u_after_rad(i) = nu_after(i) / n_after(i);
end

% ---------------------------
% Forward trace from stop to image side
% ---------------------------
y  = y_stop;
nu = nu_stop;

for i = stop_idx:N
    % State just before refraction at surface i
    y_before_mm(i) = y;
    nu_before(i)   = nu;
    u_before_rad(i)= nu / n_before(i);

    % Refraction at surface i
    nu = nu - y * phi(i);
    nu_after(i)    = nu;
    u_after_rad(i) = nu_after(i) / n_after(i);

    % Translation to the next surface (if any)
    if i < N
        y = y + (nu / n_after(i)) * t(i);
    end
end

% ---------------------------
% Object-plane position of the chief ray
% ---------------------------
% At surface 1, we know the chief ray "before" refraction: (y1, nu1)
y1  = y_before_mm(1);
nu1 = nu_before(1);

% Propagate back from surface 1 to the object plane (distance z_obj in n_air)
% y_obj_hat is the chief ray height at the object plane.
y_obj_hat = y1 - (nu1 / n_air) * z_obj;

% ---------------------------
% Pack results
% ---------------------------
res = struct();
res.Nsurf        = N;
res.R            = R;
res.t            = t;
res.n_before     = n_before;
res.n_after      = n_after;
res.phi          = phi;

res.y_before_mm  = y_before_mm;
res.nu_before    = nu_before;
res.u_before_rad = u_before_rad;
res.nu_after     = nu_after;
res.u_after_rad  = u_after_rad;

res.y_obj_hat    = y_obj_hat;   % chief-ray height at object plane
res.u_hat        = u_hat;       % chief-ray angle in object space
end
