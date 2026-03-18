import numpy as np

# const globals
gravity = 9.81  # m/s^2
kinematic_viscosity = 1.004e-6  # m^2/s for water at approx of 20 degrees Celsius
density = 998  # kg/m^3 same thing
PVC_roughness = 0.0015  # mm, convert to m in calculations

class PipeReference: # i don't know if we output from the top, center, or bottom of the pipe
    CENTERLINE = 0
    TOP = 1
    BOTTOM = -1

def effectiveHead(height, pipe_diameter, pipe_length, sin_theta, reference=PipeReference.CENTERLINE):
    offset = reference * (pipe_diameter / 2)
    return (height - offset) + pipe_length * sin_theta


def circleArea(diameter):
    return np.pi * (diameter/2)**2

def reynoldsNumber(velocity, diameter, kinematic_viscosity):
    return (velocity * diameter) / kinematic_viscosity

def frictionFactor(relative_roughness, reynolds_number, f_init = 0.02):
    # laminar
    if reynolds_number < 2000:
        return 64 / reynolds_number

    f = f_init
    for _ in range(10):
        pot_f = (-2.0 * np.log10(relative_roughness / 3.7 + 2.51 / (reynolds_number * np.sqrt(f)))) ** -2
        if abs(pot_f - f) < 1e-6:
            break
        f = pot_f

    return f

def relativeRoughness(roughness, diameter):
    return roughness / diameter

def pipeChange(height, pipe_length, sin_theta, basin_area, pipe_area, pipe_diameter, k_entry, friction_factor):
    
    coeff = pipe_area / basin_area

    denominator = 1 - (coeff ** 2) + friction_factor * (pipe_length / pipe_diameter) + k_entry
    numerator = 2 * gravity * (effectiveHead(height, pipe_diameter, pipe_length, sin_theta))

    return -coeff * np.sqrt(numerator / denominator)

def estimateFrictionFactor(height, pipe_length, sin_theta, basin_area, pipe_area, pipe_diameter, k_entry, kinematic_viscosity, f_init = 0.02):
    """
    I made this for iteratively finding the friction factor at the given height
    """
    f_guess = f_init
    coeff = pipe_area / basin_area
    new_f = 0
    while abs(f_guess - new_f) > 1e-6:
        f_guess = new_f if new_f != 0 else f_guess

        denominator = 1 - (coeff ** 2) + f_guess * (pipe_length / pipe_diameter) + k_entry
        numerator = 2 * gravity * (effectiveHead(height, pipe_diameter, pipe_length, sin_theta))
        velocity = np.sqrt(numerator / denominator)
        reynolds = reynoldsNumber(velocity, pipe_diameter, kinematic_viscosity)
        relative_roughness = relativeRoughness(PVC_roughness, pipe_diameter)
        new_f = frictionFactor(relative_roughness, reynolds)
    return new_f


def rk4(height, pipe_length, sin_theta, basin_area, pipe_area, pipe_diameter, k_entry, kinematic_viscosity, f_init=0.02, n=1000, s=0.01):
    
    def f(h):
        #friction_factor = estimateFrictionFactor(h, pipe_length, sin_theta, basin_area, pipe_area, pipe_diameter, k_entry, kinematic_viscosity, f_init)
        return pipeChange(h, pipe_length, sin_theta, basin_area, pipe_area, pipe_diameter, k_entry, 0.02)
    
    heights = [height]
    times = [0]
    curr_height = height
    curr_time = 0
    h = s  
    
    for _ in range(n):
        k1 = h * f(curr_height)
        k2 = h * f(curr_height + k1/2)
        k3 = h * f(curr_height + k2/2)
        k4 = h * f(curr_height + k3)

        curr_height += (k1 + 2*k2 + 2*k3 + k4) / 6
        curr_time += h

        heights.append(curr_height)
        times.append(curr_time)

        if height - curr_height > 0.08:
            break
    
    return times, heights


def main():
    # Basin
    initial_height = 0.1        # 10 cm from the center fo the pipe input
    basin_area = 0.32 * 0.26
    # Pipe
    pipe_diameter = 0.00794        
    pipe_length = 0.2          
    sin_theta = 1/500      

    # Loss coefficients
    k_entry = 0.5               # sharp-edged entry

    # RK4 settings
    n_steps = 10000
    step_size = 0.5            # seconds

    # Derived
    pipe_area = circleArea(pipe_diameter)

    times, heights = rk4(
        height=initial_height,
        pipe_length=pipe_length,
        sin_theta=sin_theta,
        basin_area=basin_area,
        pipe_area=pipe_area,
        pipe_diameter=pipe_diameter,
        k_entry=k_entry,
        kinematic_viscosity=kinematic_viscosity,
        n=n_steps,
        s=step_size
    )

    drop = initial_height - heights[-1]
    print(f"Simulated {times[-1]:.2f} seconds over {len(times)} steps")
    print(f"Height dropped from {initial_height:.3f} m to {heights[-1]:.3f} m ({drop*100:.1f} cm)")

    if drop >= 0.08:
        print(f"8 cm drop reached at t = {times[-1]:.2f} s")
    else:
        print(f"8 cm drop not reached in simulation window")


if __name__ == "__main__":
    main()