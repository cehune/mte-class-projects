import numpy as np

# const globals
gravity = 9.81  # m/s^2
kinematic_viscosity = 1.004e-6  # m^2/s for water at approx of 20 degrees Celsius
density = 998  # kg/m^3 same thing
PVC_roughness = 0.0015 / 1000  # mm, convert to m in calculations

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
    # laminarprint(f"Reynolds number: {reynolds_number:.2f}, Relative roughness: {relative_roughness:.6e}")
    bot = 2300
    top = 4000

    if reynolds_number < bot:
        return 64 / reynolds_number
    elif reynolds_number > top:
        return 0.25 / (np.log10(relative_roughness/3.7 + 5.74/(reynolds_number**0.9)))**2 
    print("woo");
    lam = 64 / reynolds_number
    turb = 0.25 / (np.log10(relative_roughness/3.7 + 5.74/(reynolds_number**0.9)))**2
    lerp = (reynolds_number - bot) / (top - bot)
    return lam * (1 - lerp) + turb * lerp

def relativeRoughness(roughness, diameter):
    return roughness / diameter

def pipeChange(height, pipe_length, sin_theta, basin_area, pipe_area, pipe_diameter, k_entry, k_exit, friction_factor):
    
    coeff = pipe_area / basin_area

    denominator = 1 - (coeff ** 2) + friction_factor * (pipe_length / pipe_diameter) + k_entry + k_exit
    numerator = 2 * gravity * (effectiveHead(height, pipe_diameter, pipe_length, sin_theta))

    return -coeff * np.sqrt(numerator / denominator)

def estimateFrictionFactor(height, pipe_length, sin_theta, basin_area, pipe_area, pipe_diameter, k_entry, k_exit, kinematic_viscosity, f_init=0.02):
    coeff = pipe_area / basin_area
    f = f_init

    for _ in range(50):
        L_eff = max(pipe_length - 20 * pipe_diameter, 0)
        denominator = 1 - coeff**2 + f * (pipe_length / pipe_diameter) + k_entry + k_exit
        numerator = 2 * gravity * effectiveHead(height, pipe_diameter, pipe_length, sin_theta)
        velocity = np.sqrt(numerator / denominator)
        reynolds = reynoldsNumber(velocity, pipe_diameter, kinematic_viscosity)
        rr = relativeRoughness(PVC_roughness, pipe_diameter)
        new_f = frictionFactor(rr, reynolds, f)
        #print(f"L={pipe_length:.2f} Re={reynolds:.0f} f={new_f:.4f}")

        if abs(new_f - f) < 1e-8:
            break
        f = new_f
    return f

def rk4(height, pipe_length, sin_theta, basin_area, pipe_area, pipe_diameter, k_entry, k_exit, kinematic_viscosity, f_init=0.02, n=1000, s=0.01):
    
    def f(h):
        #friction_factor = 0.04
        friction_factor = estimateFrictionFactor(h, pipe_length, sin_theta, basin_area, pipe_area, pipe_diameter, k_entry, k_exit, kinematic_viscosity, f_init)
        return pipeChange(h, pipe_length, sin_theta, basin_area, pipe_area, pipe_diameter, k_entry, k_exit, friction_factor)
    
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
    initial_height = 0.1        # 10 cm from the center fo the pipe input
    basin_area = 0.32 * 0.26
    # Pipe
    pipe_diameter = 0.00794        
    sin_theta = 1/150      

    # Loss coefficients
    k_entry = 0.6               # sharp-edged entry
    k_exit = 0.0  # claude said i'm double counting with the 1 in my denom already lol need to verify
    # RK4 settings
    n_steps = 10000
    step_size = 0.1            # seconds

    # Derived
    pipe_area = circleArea(pipe_diameter)

    for pipe_length in [0.2, 0.3, 0.4, 0.6]:
        times, heights = rk4(
            height=initial_height,
            pipe_length=pipe_length,
            sin_theta=sin_theta,
            basin_area=basin_area,
            pipe_area=pipe_area,
            pipe_diameter=pipe_diameter,
            k_entry=k_entry,
            k_exit=k_exit,
            kinematic_viscosity=kinematic_viscosity,
            n=n_steps,
            s=step_size
        )

        drop = initial_height - heights[-1]
        print(f"\nPipe length: {pipe_length:.2f} m")
        # print(f"Simulated {times[-1]:.2f} seconds over {len(times)} steps")
        # print(f"Height dropped from {initial_height:.3f} m to {heights[-1]:.3f} m ({drop*100:.1f} cm)")

        if drop >= 0.08:
            print(f"8 cm drop reached at t = {times[-1]:.2f} s\n")
        #else:
            #print(f"8 cm drop not reached in simulation window\n")
        


if __name__ == "__main__":
    main()