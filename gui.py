import tkinter as tk
from tkinter import ttk

# For embedding matplotlib in Tkinter:
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# We assume you have "foc_sim.py" in the same folder with:
#   default_params, run_foc_simulation
# containing all the motor, load, and control logic.
from foc_sim import default_params, run_foc_simulation

###############################################################################
# A helper function for creating a tooltip on any widget
###############################################################################
def create_tooltip(widget, text):
    """
    Attach a tooltip to the given widget. 
    'text' is the string displayed in the tooltip popup.
    """
    def on_enter(event):
        # Create a Toplevel window (no window decorations)
        widget.tooltip = tk.Toplevel(widget)
        widget.tooltip.overrideredirect(True)  # no title bar, borders, etc.

        # Position it near the widget
        x = widget.winfo_rootx() + 50
        y = widget.winfo_rooty() + 20
        widget.tooltip.geometry(f"+{x}+{y}")

        # Create a label in that Toplevel
        label = tk.Label(widget.tooltip, text=text,
                         background="lightyellow", relief="solid",
                         borderwidth=1, padx=5, pady=2)
        label.pack()

    def on_leave(event):
        if hasattr(widget, 'tooltip'):
            widget.tooltip.destroy()
            widget.tooltip = None

    # Bind mouse events
    widget.bind("<Enter>", on_enter)
    widget.bind("<Leave>", on_leave)

###############################################################################
# The main GUI class
###############################################################################
class FOCGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FOC Simulation GUI with Combined Parameters Tab + Tooltips")

        # We'll keep a copy of the user parameters
        self.user_params = dict(default_params)

        # We'll store trajectory settings
        self.angle_final_deg = tk.StringVar(value="30.0")
        self.v_max = tk.StringVar(value="10.0")
        self.a_max = tk.StringVar(value="20.0")
        self.gear_ratio = tk.StringVar(value="3.0")

        # Create a Notebook with 2 tabs: "Parameters" and "Plot"
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # --- Tab 1: Parameters ---
        self.params_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.params_tab, text="Parameters")

        # Inside this tab, we'll place 4 labeled frames side-by-side (columns):
        #  1) Motor & Load
        #  2) Controller Gains
        #  3) Simulation
        #  4) Trajectory

        self.params_tab.columnconfigure(0, weight=1)
        self.params_tab.columnconfigure(1, weight=1)
        self.params_tab.columnconfigure(2, weight=1)
        self.params_tab.columnconfigure(3, weight=1)

        # Build each column
        self._build_motor_frame(self.params_tab, 0)  # col 0
        self._build_gains_frame(self.params_tab, 1)  # col 1
        self._build_sim_frame(self.params_tab, 2)    # col 2
        self._build_traj_frame(self.params_tab, 3)   # col 3

        # --- Tab 2: Plot ---
        self.plot_tab = ttk.Frame(self.notebook)
        self.notebook.add(self.plot_tab, text="Plot")

        # Create a matplotlib Figure + Canvas in the "Plot" tab
        self.figure = Figure(figsize=(8,6))
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.plot_tab)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Button to run the simulation
        run_button = ttk.Button(self, text="Run Simulation", command=self.run_simulation)
        run_button.pack(pady=5)

    # -------------------------------------------------------------------------
    # 1) MOTOR & LOAD FRAME
    # -------------------------------------------------------------------------
    def _build_motor_frame(self, parent, col_index):
        """Create a LabelFrame in the given column for Motor & Load parameters."""
        lf = ttk.LabelFrame(parent, text="Motor & Load")
        lf.grid(row=0, column=col_index, padx=5, pady=5, sticky="nsew")

        row_idx = 0
        # R_s
        label_Rs = ttk.Label(lf, text="R_s [ohms]")
        label_Rs.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_Rs, 
            "Stator phase resistance (per phase).\n"
            "e.g., measure with a multimeter or from datasheet.")
        
        self.entry_Rs = ttk.Entry(lf, width=10)
        self.entry_Rs.insert(0, str(self.user_params['R_s']))
        self.entry_Rs.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # L
        label_L = ttk.Label(lf, text="L [H]")
        label_L.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_L, 
            "Stator inductance (assuming Ld = Lq) in Henries.\n"
            "From datasheet or low-frequency measurement.")
        
        self.entry_L = ttk.Entry(lf, width=10)
        self.entry_L.insert(0, str(self.user_params['L']))
        self.entry_L.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # p
        label_p = ttk.Label(lf, text="p (pole pairs)")
        label_p.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_p, 
            "Number of pole pairs, e.g. if the motor\n"
            "has 8 total poles, then p=4.")
        
        self.entry_p = ttk.Entry(lf, width=10)
        self.entry_p.insert(0, str(self.user_params['p']))
        self.entry_p.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # lambda_m
        label_lam = ttk.Label(lf, text="lambda_m [Wb]")
        label_lam.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_lam, 
            "PM flux linkage. Often derived from the\n"
            "back-EMF constant. e.g. 0.015 Wb.")
        
        self.entry_lam = ttk.Entry(lf, width=10)
        self.entry_lam.insert(0, str(self.user_params['lambda_m']))
        self.entry_lam.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # J_m
        label_Jm = ttk.Label(lf, text="J_m [kg*m^2]")
        label_Jm.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_Jm, 
            "Rotor (motor) inertia.\n"
            "From datasheet or spin-down tests.")
        
        self.entry_Jm = ttk.Entry(lf, width=10)
        self.entry_Jm.insert(0, str(self.user_params['J_m']))
        self.entry_Jm.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # B_m
        label_Bm = ttk.Label(lf, text="B_m [N*m*s/rad]")
        label_Bm.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_Bm, 
            "Viscous friction for the motor side.\n"
            "Might be small if the motor is well-lubricated.")
        
        self.entry_Bm = ttk.Entry(lf, width=10)
        self.entry_Bm.insert(0, str(self.user_params['B_m']))
        self.entry_Bm.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # J_load
        label_Jload = ttk.Label(lf, text="J_load [kg*m^2]")
        label_Jload.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_Jload, 
            "Load inertia at the load shaft.\n"
            "Calculate from geometry or measure.")
        
        self.entry_Jload = ttk.Entry(lf, width=10)
        self.entry_Jload.insert(0, str(self.user_params['J_load']))
        self.entry_Jload.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # B_load
        label_Bload = ttk.Label(lf, text="B_load [N*m*s/rad]")
        label_Bload.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_Bload,
            "Viscous friction for the load side.\n"
            "If unknown, start with a small guess or measure.")
        
        self.entry_Bload = ttk.Entry(lf, width=10)
        self.entry_Bload.insert(0, str(self.user_params['B_load']))
        self.entry_Bload.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # T_load
        label_Tload = ttk.Label(lf, text="T_load [Nm]")
        label_Tload.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_Tload, 
            "Constant external load torque.\n"
            "If gravity or friction exist, put it here.")
        
        self.entry_Tload = ttk.Entry(lf, width=10)
        self.entry_Tload.insert(0, str(self.user_params['T_load']))
        self.entry_Tload.grid(row=row_idx, column=1, padx=5, pady=2)

    # -------------------------------------------------------------------------
    # 2) CONTROLLER GAINS FRAME
    # -------------------------------------------------------------------------
    def _build_gains_frame(self, parent, col_index):
        lf = ttk.LabelFrame(parent, text="Controller Gains")
        lf.grid(row=0, column=col_index, padx=5, pady=5, sticky="nsew")

        row_idx = 0

        # Kp_pos
        label_Kp_pos = ttk.Label(lf, text="Kp_pos")
        label_Kp_pos.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_Kp_pos,
            "Proportional gain for the position loop.\n"
            "Higher Kp_pos => faster response, but can overshoot.")
        
        self.entry_Kp_pos = ttk.Entry(lf, width=8)
        self.entry_Kp_pos.insert(0, str(self.user_params['Kp_pos']))
        self.entry_Kp_pos.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # Ki_pos
        label_Ki_pos = ttk.Label(lf, text="Ki_pos")
        label_Ki_pos.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_Ki_pos,
            "Integral gain for the position loop.\n"
            "Helps remove steady-state error but can cause overshoot if too high.")
        
        self.entry_Ki_pos = ttk.Entry(lf, width=8)
        self.entry_Ki_pos.insert(0, str(self.user_params['Ki_pos']))
        self.entry_Ki_pos.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # Kd_pos
        label_Kd_pos = ttk.Label(lf, text="Kd_pos")
        label_Kd_pos.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_Kd_pos,
            "Derivative gain for position.\n"
            "Small value for damping overshoot.")
        
        self.entry_Kd_pos = ttk.Entry(lf, width=8)
        self.entry_Kd_pos.insert(0, str(self.user_params['Kd_pos']))
        self.entry_Kd_pos.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # Kp_speed
        label_Kp_speed = ttk.Label(lf, text="Kp_speed")
        label_Kp_speed.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_Kp_speed,
            "Proportional gain for the speed loop.\n"
            "If too high => potential oscillations.")
        
        self.entry_Kp_speed = ttk.Entry(lf, width=8)
        self.entry_Kp_speed.insert(0, str(self.user_params['Kp_speed']))
        self.entry_Kp_speed.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # Ki_speed
        label_Ki_speed = ttk.Label(lf, text="Ki_speed")
        label_Ki_speed.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_Ki_speed,
            "Integral gain for speed loop.\n"
            "Boosts torque at steady-state error, watch out for windup.")
        
        self.entry_Ki_speed = ttk.Entry(lf, width=8)
        self.entry_Ki_speed.insert(0, str(self.user_params['Ki_speed']))
        self.entry_Ki_speed.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # Kp_id
        label_Kp_id = ttk.Label(lf, text="Kp_id")
        label_Kp_id.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_Kp_id,
            "Proportional gain for the d-axis current loop.\n"
            "Typically set to get fast current response.")
        
        self.entry_Kp_id = ttk.Entry(lf, width=8)
        self.entry_Kp_id.insert(0, str(self.user_params['Kp_id']))
        self.entry_Kp_id.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # Ki_id
        label_Ki_id = ttk.Label(lf, text="Ki_id")
        label_Ki_id.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_Ki_id,
            "Integral gain for d-axis current loop.\n"
            "Counteracts steady-state errors in d-axis current.")
        
        self.entry_Ki_id = ttk.Entry(lf, width=8)
        self.entry_Ki_id.insert(0, str(self.user_params['Ki_id']))
        self.entry_Ki_id.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # Kp_iq
        label_Kp_iq = ttk.Label(lf, text="Kp_iq")
        label_Kp_iq.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_Kp_iq,
            "Proportional gain for the q-axis current loop.\n"
            "Mainly controls torque response.")
        
        self.entry_Kp_iq = ttk.Entry(lf, width=8)
        self.entry_Kp_iq.insert(0, str(self.user_params['Kp_iq']))
        self.entry_Kp_iq.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # Ki_iq
        label_Ki_iq = ttk.Label(lf, text="Ki_iq")
        label_Ki_iq.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_Ki_iq,
            "Integral gain for q-axis current loop.\n"
            "Helps maintain correct torque under load variations.")
        
        self.entry_Ki_iq = ttk.Entry(lf, width=8)
        self.entry_Ki_iq.insert(0, str(self.user_params['Ki_iq']))
        self.entry_Ki_iq.grid(row=row_idx, column=1, padx=5, pady=2)

    # -------------------------------------------------------------------------
    # 3) SIMULATION FRAME
    # -------------------------------------------------------------------------
    def _build_sim_frame(self, parent, col_index):
        lf = ttk.LabelFrame(parent, text="Simulation")
        lf.grid(row=0, column=col_index, padx=5, pady=5, sticky="nsew")

        row_idx = 0

        # t_sim
        label_tsim = ttk.Label(lf, text="t_sim [s]")
        label_tsim.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_tsim,
            "Total simulation time in seconds.\n"
            "Make sure it's long enough for your motion to complete.")
        
        self.entry_tsim = ttk.Entry(lf, width=8)
        self.entry_tsim.insert(0, str(self.user_params['t_sim']))
        self.entry_tsim.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # dt
        label_dt = ttk.Label(lf, text="dt [s]")
        label_dt.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_dt,
            "Simulation time-step (seconds).\n"
            "Smaller => more accurate, but slower simulation.")
        
        self.entry_dt = ttk.Entry(lf, width=8)
        self.entry_dt.insert(0, str(self.user_params['dt']))
        self.entry_dt.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # V_dc
        label_vdc = ttk.Label(lf, text="V_dc [V]")
        label_vdc.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_vdc,
            "DC bus (supply) voltage for your inverter.\n"
            "e.g. 12 V, 24 V, or higher.")
        
        self.entry_vdc = ttk.Entry(lf, width=8)
        self.entry_vdc.insert(0, str(self.user_params['V_dc']))
        self.entry_vdc.grid(row=row_idx, column=1, padx=5, pady=2)

    # -------------------------------------------------------------------------
    # 4) TRAJECTORY FRAME
    # -------------------------------------------------------------------------
    def _build_traj_frame(self, parent, col_index):
        lf = ttk.LabelFrame(parent, text="Trajectory")
        lf.grid(row=0, column=col_index, padx=5, pady=5, sticky="nsew")

        row_idx = 0

        # Final Angle
        label_angle = ttk.Label(lf, text="Final Angle (deg)")
        label_angle.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_angle,
            "The end position for the load in degrees.\n"
            "Converted to radians internally.")
        
        entry_angle = ttk.Entry(lf, textvariable=self.angle_final_deg, width=8)
        entry_angle.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # v_max
        label_vmax = ttk.Label(lf, text="v_max (rad/s)")
        label_vmax.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_vmax,
            "Maximum velocity for the trapezoidal profile.\n"
            "If distance is short, might not be reached.")
        
        entry_vmax = ttk.Entry(lf, textvariable=self.v_max, width=8)
        entry_vmax.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # a_max
        label_amax = ttk.Label(lf, text="a_max (rad/s^2)")
        label_amax.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_amax,
            "Maximum acceleration for the trapezoidal profile.\n"
            "Higher => faster moves, but more torque required.")
        
        entry_amax = ttk.Entry(lf, textvariable=self.a_max, width=8)
        entry_amax.grid(row=row_idx, column=1, padx=5, pady=2)
        row_idx += 1

        # gear_ratio
        label_gear = ttk.Label(lf, text="Gear Ratio")
        label_gear.grid(row=row_idx, column=0, sticky=tk.W)
        create_tooltip(label_gear,
            "Ratio between motor and load speeds.\n"
            "motor_speed = gear_ratio * load_speed.")
        
        entry_gear = ttk.Entry(lf, textvariable=self.gear_ratio, width=8)
        entry_gear.grid(row=row_idx, column=1, padx=5, pady=2)

    # -------------------------------------------------------------------------
    # RUN SIM
    # -------------------------------------------------------------------------
    def run_simulation(self):
        """Read inputs, run simulation, then plot in the 'Plot' tab."""
        try:
            # 1) Gather from Motor & Load
            self.user_params['R_s']      = float(self.entry_Rs.get())
            self.user_params['L']        = float(self.entry_L.get())
            self.user_params['p']        = float(self.entry_p.get())
            self.user_params['lambda_m'] = float(self.entry_lam.get())
            self.user_params['J_m']      = float(self.entry_Jm.get())
            self.user_params['B_m']      = float(self.entry_Bm.get())
            self.user_params['J_load']   = float(self.entry_Jload.get())
            self.user_params['B_load']   = float(self.entry_Bload.get())
            self.user_params['T_load']   = float(self.entry_Tload.get())

            # 2) Gains
            self.user_params['Kp_pos']    = float(self.entry_Kp_pos.get())
            self.user_params['Ki_pos']    = float(self.entry_Ki_pos.get())
            self.user_params['Kd_pos']    = float(self.entry_Kd_pos.get())
            self.user_params['Kp_speed']  = float(self.entry_Kp_speed.get())
            self.user_params['Ki_speed']  = float(self.entry_Ki_speed.get())
            self.user_params['Kp_id']     = float(self.entry_Kp_id.get())
            self.user_params['Ki_id']     = float(self.entry_Ki_id.get())
            self.user_params['Kp_iq']     = float(self.entry_Kp_iq.get())
            self.user_params['Ki_iq']     = float(self.entry_Ki_iq.get())

            # 3) Simulation
            self.user_params['t_sim']     = float(self.entry_tsim.get())
            self.user_params['dt']        = float(self.entry_dt.get())
            self.user_params['V_dc']      = float(self.entry_vdc.get())

            # 4) Trajectory
            angle_deg  = float(self.angle_final_deg.get())
            v_max      = float(self.v_max.get())
            a_max      = float(self.a_max.get())
            gear_ratio = float(self.gear_ratio.get())

        except ValueError:
            print("Invalid input. Please check fields.")
            return

        # Run the simulation
        results = run_foc_simulation(
            params=self.user_params,
            angle_final_deg=angle_deg,
            v_max=v_max,
            a_max=a_max,
            gear_ratio=gear_ratio
        )

        # Plot the results inside the 'Plot' tab
        self.draw_plots_in_figure(results)

        # Switch to the "Plot" tab
        plot_tab_idx = self.notebook.index(self.plot_tab)
        self.notebook.select(plot_tab_idx)

    def draw_plots_in_figure(self, results):
        """Draw the 4x2 subplots in the existing self.figure."""
        self.figure.clear()
        axs = self.figure.subplots(4, 2, sharex=True)
        axs = axs.ravel()
        t = results['time']

        # Example: Set human reaction time (in seconds)
        reaction_time = 0.25  # you can choose any value

        # 0) Load Angle vs. Reference
        axs[0].plot(t, results['theta_load'], label='Load Angle (rad)')
        axs[0].plot(t, results['ref_angle'], 'r--', label='Ref Angle (rad)')
        axs[0].axvline(reaction_time, color='limegreen', linestyle='--', label='Human Reaction')
        axs[0].set_title("Load Angle vs. Reference")
        axs[0].legend()
        axs[0].grid(True)

        # 1) Speeds
        axs[1].plot(t, results['omega_load'], label='Load Speed')
        axs[1].plot(t, results['omega_motor'], label='Motor Speed')
        axs[1].set_title("Speeds (rad/s)")
        axs[1].legend()
        axs[1].grid(True)

        # 2) dq Currents
        axs[2].plot(t, results['i_d'], label='i_d')
        axs[2].plot(t, results['i_q'], label='i_q')
        axs[2].set_title("dq Currents (A)")
        axs[2].legend()
        axs[2].grid(True)

        # 3) dq Voltages
        axs[3].plot(t, results['v_d'], label='v_d')
        axs[3].plot(t, results['v_q'], label='v_q')
        axs[3].set_title("dq Voltages (V)")
        axs[3].legend()
        axs[3].grid(True)

        # 4) Electromagnetic Torque
        axs[4].plot(t, results['torque_e'], label='Torque (Nm)')
        axs[4].set_title("Electromagnetic Torque")
        axs[4].legend()
        axs[4].grid(True)

        # 5) Phase Voltages
        axs[5].plot(t, results['v_a'], label='v_a')
        axs[5].plot(t, results['v_b'], label='v_b')
        axs[5].plot(t, results['v_c'], label='v_c')
        axs[5].set_title("Three-Phase Voltages")
        axs[5].legend()
        axs[5].grid(True)

        # 6) Duty Cycles
        axs[6].plot(t, results['d_a'], label='d_a')
        axs[6].plot(t, results['d_b'], label='d_b')
        axs[6].plot(t, results['d_c'], label='d_c')
        axs[6].set_title("SVPWM Duty Cycles")
        axs[6].legend()
        axs[6].grid(True)

        # 7) Spare
        axs[7].text(0.5, 0.5, "Spare Plot\n(Adjust as needed)",
                    ha='center', va='center', fontsize=12,
                    bbox=dict(facecolor='lightgray', alpha=0.5))
        axs[7].set_axis_off()

        for ax in axs:
            ax.set_xlabel("Time (s)")

        self.figure.tight_layout()
        self.canvas.draw()

def main():
    app = FOCGUI()
    app.mainloop()

if __name__ == "__main__":
    main()
