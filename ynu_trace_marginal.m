function res = ynu_trace_marginal(surfaces, n_obj, z_to_surf1, u0)
% Computes the marginal-ray Y–ν trace through a sequence of refracting
% surfaces using paraxial (first-order) optics.
%
% surfaces: [R, t_to_next, n_after]  — an N×3 description of the lens system
% n_obj: refractive index before the first surface
% z_to_surf1: axial distance from object plane to surface 1, in index n_obj
% u0: initial paraxial angle in object space (rad)

N = size(surfaces, 1);
R = surfaces(:, 1);   % surface radii
t = surfaces(:, 2);   % axial spacings to next surface
n_after = surfaces(:, 3);   % refractive index after each surface

% Allocate arrays for ray table
y_before = zeros(N, 1);   % height before refraction
u_before = zeros(N, 1);   % slope before refraction
nu_after = zeros(N, 1);   % optical direction cosine after refraction
n_before = zeros(N, 1);   % refractive index before surface
y_after  = zeros(N, 1);   % height after refraction
phi_list = zeros(N, 1);   % surface powers

% Initial conditions in object space
n_curr = n_obj;
y = u0 * z_to_surf1;   % marginal ray height at surface 1
u = u0;

for i = 1:N
    % State just before refraction at surface i
    y_before(i) = y;
    u_before(i) = u;
    n_before(i) = n_curr;
    nu_before = n_curr * u;

    % Surface power phi
    if isinf(R(i))
        phi = 0;   % plane or stop -> no refraction power
    else
        phi = (n_after(i) - n_curr) / R(i);
    end
    phi_list(i) = phi;

    % Refraction: nu_after = nu_before − y*phi
    nu_a = nu_before - y * phi;
    nu_after(i) = nu_a;

    % Propagate through thickness to next surface
    n_next = n_after(i);
    y_after(i) = y;
    y = y + (nu_a / n_next) * t(i);
    u = nu_a / n_next;

    % Prepare index for next interval
    n_curr = n_next;
end

% ---- Image plane calculation after final surface ----

n_last = n_curr;
y_last = y_before(end);
nu_last = nu_after(end);

% Paraxial image distance from last surface:  s' = -n·y / ν
s_img = - n_last * y_last / nu_last;

% Collect results
res = struct();
res.Nsurf      = N;
res.R          = R;
res.t          = t;
res.n_before   = n_before;
res.n_after    = n_after;
res.y_before   = y_before;
res.y_after    = y_after;
res.u_before   = u_before;
res.nu_after   = nu_after;
res.phi        = phi_list;

res.n_last     = n_last;
res.y_last     = y_last;
res.nu_last    = nu_last;
res.s_img      = s_img;
end
