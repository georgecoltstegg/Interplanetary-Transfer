% given
mu_s = 130000000000; %km^3/s^2
mu_L = 355000; % km^3/s^2
AU = 149597870.691; %km
R_A = 1*AU; %km
R_L = 0.75*AU; %km
r_L = 5500; %km
theta_vec = linspace(-39, -30, 10); % deg
Z_p_vec = linspace(300, 400, 51); % km
dV_matrix = zeros(length(theta_vec));
iterations = 0;

for i = 1:length(theta_vec)
    % Find V1v pre fly-by ellipse
    theta = theta_vec(i);
    e_1 = (R_A - R_L)/(R_A + R_L*cos(theta*(pi/180)));
    h_1 = sqrt(R_A*mu_s*(1-e_1));
    v_tan = h_1/R_L; %km/s
    v_rad = (mu_s/h_1)*e_1*sin(theta*(pi/180)); %km/s
    V1v_vector = [v_tan -v_rad]; % putting rad and tan components together
    V1v = sqrt(v_tan^2 + v_rad^2); %km/s magnitude of V1v vector
    V_cL = sqrt(mu_s/R_L); %circular orbit velocity of planet longhorn
    V_cL_vector = [V_cL 0]; % circular orbit in vector form
    v_inf1_vector = V1v_vector - V_cL_vector;
    v_inf = sqrt(dot(v_inf1_vector, v_inf1_vector)); %finding v_inf
    % using these values for theta, iterate through every altitude
    for j = 1:length(Z_p_vec)
        Z_p = Z_p_vec(j); % altitude
        r_p = r_L + Z_p; %radius
        h = r_p*sqrt(v_inf^2 + ((2*mu_L)/r_p));
        e = 1 + ((r_p*v_inf^2)/mu_L); %hyperbolic eccentricity
        turn = 2*asin(1/e); % turn angle
        theta_inf = acos(-(1/e)); %hyperbolic true anomaly
        phi_1 = atan(v_inf1_vector(2) / v_inf1_vector(1));
        phi_2 = phi_1 - turn;
        v_inf2__vector = v_inf*[cos(phi_2) sin(phi_2)];
        V2v_vector = V_cL_vector + v_inf2__vector;
        V2v = sqrt((V2v_vector(1))^2 + (V2v_vector(2))^2); % second impulse
        % calculate delta V at specific theta and altitude
        delta_V = abs(V2v - V1v);
        % add delta v to correct position in matrix
        dV_matrix(i, j) = delta_V;
        % count total iterations
        iterations = iterations + 1;
    end
end

% plot contour map
figure
contourf(Z_p_vec, theta_vec, dV_matrix)
colormap jet
colorbar('vertical');
ylabel('True Anomaly [deg]');xlabel('Altitude [km]')
title('Impulse Contour [km/s]');set(gca, 'fontsize', 18)
