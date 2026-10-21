%Rate function describing Newton's law of gravitation
%INPUTS:
%t: the current time
%V: the vector of the position and velocity of the planet
% V = [x_p; y_p; dxdt_p; dydt_p]
%orbit_params: a struct describing the system parameters
% orbit_params.m_sun is the mass of the sun
% orbit_params.m_planet is the mass of the planet
% orbit_params.G is the gravitational constant
%OUTPUTS:
%dVdt: a column vector describing the time derivative of V:
% dVdt = [dxdt_p; dydt_p; d2xdt2_p; d2ydt2_p]
function dVdt = gravity_rate_func(t,V,orbit_params)
    pos_x = V(1);
    pos_y = V(2);
    vel_x = V(3);
    vel_y = V(4);

    ms = orbit_params.m_sun;
    G = orbit_params.G;

    r = sqrt(pos_x^2 + pos_y^2);

    accel_x = -G*ms*pos_x/r^3;
    accel_y = -G*ms*pos_y/r^3;

    dVdt = [vel_x; vel_y; accel_x; accel_y];
end