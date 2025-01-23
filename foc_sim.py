# foc_sim.py
import numpy as np

default_params = {
    'R_s': 0.5,
    'L': 1e-3,
    'p': 4,
    'lambda_m': 0.015,
    'J_m': 1e-4,
    'B_m': 1e-4,
    'J_load': 5e-4,
    'B_load': 2e-4,
    'T_load': 0.0,
    'Kp_pos': 2.0,
    'Ki_pos': 0.3,
    'Kd_pos': 0.05,
    'Kp_speed': 0.02,
    'Ki_speed': 0.2,
    'Kp_id': 20.0,
    'Ki_id': 2000.0,
    'Kp_iq': 20.0,
    'Ki_iq': 2000.0,
    't_sim': 2.0,
    'dt': 1e-5,
    'V_dc': 24.0,
}

def angle_trapezoid_profile(t, distance, v_max, a_max):
    # (same trapezoid code)
    sign = 1.0
    if distance < 0:
        sign = -1.0
        distance = abs(distance)
    t_accel = v_max / a_max
    dist_accel = 0.5 * a_max * (t_accel**2)
    if 2.0 * dist_accel > distance:
        # Triangular
        import math
        T = np.sqrt(2.0 * distance / a_max)
        half_T = 0.5 * T
        if t < 0:
            pos = 0
        elif t < half_T:
            pos = 0.5*a_max*t**2
        elif t < T:
            td = t - half_T
            dist_first = 0.5*a_max*(half_T**2)
            pos = dist_first + (a_max*half_T*td - 0.5*a_max*td**2)
        else:
            pos = distance
    else:
        # Full trapezoid
        dist_cruise = distance - 2.0*dist_accel
        t_cruise = dist_cruise / v_max
        t_total = t_accel + t_cruise + t_accel
        if t < 0:
            pos = 0.0
        elif t < t_accel:
            pos = 0.5 * a_max * t**2
        elif t < (t_accel + t_cruise):
            dist_first = dist_accel
            pos = dist_first + v_max*(t - t_accel)
        elif t < t_total:
            td = t - (t_accel + t_cruise)
            dist_first = dist_accel + dist_cruise
            pos = dist_first + (v_max*td - 0.5*a_max*td**2)
        else:
            pos = distance
    return sign * pos


def compute_svpwm_duties(v_d, v_q, theta_e, V_dc):
    # (same inverse Park + Clarke code, no plt)
    import numpy as np
    cos_e = np.cos(theta_e)
    sin_e = np.sin(theta_e)

    v_alpha = cos_e*v_d - sin_e*v_q
    v_beta  = sin_e*v_d + cos_e*v_q
    v_a = v_alpha
    v_b = -0.5*v_alpha + (np.sqrt(3)/2)*v_beta
    v_c = -0.5*v_alpha - (np.sqrt(3)/2)*v_beta

    vs = np.array([v_a, v_b, v_c])
    v_min = vs.min()
    v_max = vs.max()

    offset = -0.5*(v_max + v_min)
    v_a_shift = v_a + offset
    v_b_shift = v_b + offset
    v_c_shift = v_c + offset

    d_a = v_a_shift / V_dc
    d_b = v_b_shift / V_dc
    d_c = v_c_shift / V_dc

    d_a = np.clip(d_a, 0.0, 1.0)
    d_b = np.clip(d_b, 0.0, 1.0)
    d_c = np.clip(d_c, 0.0, 1.0)

    return (v_a, v_b, v_c), (d_a, d_b, d_c)


def run_foc_simulation(params, angle_final_deg=30.0, v_max=10.0, a_max=20.0, gear_ratio=3.0):
    import numpy as np
    import math

    R_s       = params['R_s']
    L         = params['L']
    p         = params['p']
    lambda_m  = params['lambda_m']
    J_m       = params['J_m']
    B_m       = params['B_m']
    J_load    = params['J_load']
    B_load    = params['B_load']
    T_load    = params['T_load']

    Kp_pos    = params['Kp_pos']
    Ki_pos    = params['Ki_pos']
    Kd_pos    = params['Kd_pos']
    Kp_speed  = params['Kp_speed']
    Ki_speed  = params['Ki_speed']
    Kp_id     = params['Kp_id']
    Ki_id     = params['Ki_id']
    Kp_iq     = params['Kp_iq']
    Ki_iq     = params['Ki_iq']

    t_sim     = params['t_sim']
    dt        = params['dt']
    n_steps   = int(t_sim/dt)
    V_dc      = params['V_dc']

    distance_rad = math.radians(angle_final_deg)

    time_arr       = np.zeros(n_steps)
    theta_load_arr = np.zeros(n_steps)
    omega_load_arr = np.zeros(n_steps)
    theta_m_arr    = np.zeros(n_steps)
    omega_m_arr    = np.zeros(n_steps)
    i_d_arr        = np.zeros(n_steps)
    i_q_arr        = np.zeros(n_steps)
    v_d_arr        = np.zeros(n_steps)
    v_q_arr        = np.zeros(n_steps)
    torque_m_arr   = np.zeros(n_steps)

    v_a_arr        = np.zeros(n_steps)
    v_b_arr        = np.zeros(n_steps)
    v_c_arr        = np.zeros(n_steps)
    d_a_arr        = np.zeros(n_steps)
    d_b_arr        = np.zeros(n_steps)
    d_c_arr        = np.zeros(n_steps)
    ref_angle_arr  = np.zeros(n_steps)

    # Integrators
    pos_int_err   = 0.0
    pos_err_prev  = 0.0
    speed_int_err = 0.0
    id_int_err    = 0.0
    iq_int_err    = 0.0

    # IC
    theta_load = 0.0
    omega_load = 0.0
    theta_m    = 0.0
    omega_m    = 0.0
    i_d        = 0.0
    i_q        = 0.0

    for k in range(n_steps):
        t = k * dt
        time_arr[k] = t

        # trapezoid
        angle_ref = angle_trapezoid_profile(t, distance_rad, v_max, a_max)
        ref_angle_arr[k] = angle_ref

        # position PID
        pos_err = angle_ref - theta_load
        pos_int_err += pos_err * dt
        pos_der_err = (pos_err - pos_err_prev)/dt
        pos_err_prev = pos_err

        omega_load_ref = (Kp_pos*pos_err + Ki_pos*pos_int_err + Kd_pos*pos_der_err)

        # speed loop
        speed_err = (gear_ratio * omega_load_ref) - omega_m
        speed_int_err += speed_err * dt
        iq_ref = (Kp_speed*speed_err + Ki_speed*speed_int_err)
        id_ref = 0.0

        # current loops
        id_err = id_ref - i_d
        iq_err = iq_ref - i_q
        id_int_err += id_err * dt
        iq_int_err += iq_err * dt

        v_d_ref = Kp_id*id_err + Ki_id*id_int_err
        v_q_ref = Kp_iq*iq_err + Ki_iq*iq_int_err

        v_d_arr[k] = v_d_ref
        v_q_arr[k] = v_q_ref

        # electrical
        omega_e = p * omega_m
        di_d_dt = (v_d_ref + omega_e*L*i_q - R_s*i_d) / L
        di_q_dt = (v_q_ref - omega_e*L*i_d - omega_e*lambda_m - R_s*i_q) / L

        i_d += di_d_dt*dt
        i_q += di_q_dt*dt

        # mechanical
        T_e = 1.5 * p * lambda_m * i_q
        domega_m_dt = (T_e - B_m*omega_m)/J_m
        omega_m += domega_m_dt*dt
        theta_m += omega_m*dt

        T_out = gear_ratio * T_e
        domega_load_dt = (T_out - B_load*omega_load - T_load)/J_load
        omega_load += domega_load_dt*dt
        theta_load += omega_load*dt

        # log
        theta_load_arr[k] = theta_load
        omega_load_arr[k] = omega_load
        theta_m_arr[k]    = theta_m
        omega_m_arr[k]    = omega_m
        i_d_arr[k]        = i_d
        i_q_arr[k]        = i_q
        torque_m_arr[k]   = T_e

        theta_e = p * theta_m
        (v_a, v_b, v_c), (d_a, d_b, d_c) = compute_svpwm_duties(v_d_ref, v_q_ref, theta_e, V_dc)
        v_a_arr[k] = v_a
        v_b_arr[k] = v_b
        v_c_arr[k] = v_c
        d_a_arr[k] = d_a
        d_b_arr[k] = d_b
        d_c_arr[k] = d_c

    results = {
        'time':         time_arr,
        'theta_load':   theta_load_arr,
        'omega_load':   omega_load_arr,
        'theta_motor':  theta_m_arr,
        'omega_motor':  omega_m_arr,
        'i_d':          i_d_arr,
        'i_q':          i_q_arr,
        'v_d':          v_d_arr,
        'v_q':          v_q_arr,
        'torque_e':     torque_m_arr,
        'v_a':          v_a_arr,
        'v_b':          v_b_arr,
        'v_c':          v_c_arr,
        'd_a':          d_a_arr,
        'd_b':          d_b_arr,
        'd_c':          d_c_arr,
        'ref_angle':    ref_angle_arr,
    }
    return results
