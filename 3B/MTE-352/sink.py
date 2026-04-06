import numpy as np

# const globals
gravity = 9.81  # m/s^2
kinematic_viscosity = 1.004e-6  # m^2/s for water at approx of 20 degrees Celsius
density = 998  # kg/m^3 same thing
PVC_roughness = 0.0015 / 1000  # mm, convert to m

class PipeReference:
    CENTERLINE = 0
    TOP = 1
    BOTTOM = -1

class Config:
    USE_T_JUNCTION = False
    K_junction = 2.0      # tune: straight-run ~0.3, branch turn ~1.0
    k_entry = 0.5         # tune: well-rounded ~0.1, sharp-edged ~0.5, re-entrant ~0.8
    initial_height = 0.1
    basin_area = 0.32 * 0.26
    pipe_diameter = 0.00794
    sin_theta = 1 / 150
    n_steps = 10000
    step_size = 0.1
    drop_target = 0.08

def effectiveHead(height, pipe_diameter, pipe_length, sin_theta, reference=PipeReference.CENTERLINE):
    offset = reference * (pipe_diameter / 2)
    return (height - offset) + pipe_length * sin_theta

def circleArea(diameter):
    return np.pi * (diameter / 2) ** 2

def reynoldsNumber(velocity, diameter, kinematic_viscosity):
    return (velocity * diameter) / kinematic_viscosity

def frictionFactor(relative_roughness, reynolds_number):
    bot = 2300
    top = 4000

    if reynolds_number < bot:
        return 64 / reynolds_number
    elif reynolds_number > top:
        return 0.25 / (np.log10(relative_roughness / 3.7 + 5.74 / (reynolds_number ** 0.9))) ** 2

    lam = 64 / reynolds_number
    turb = 0.25 / (np.log10(relative_roughness / 3.7 + 5.74 / (reynolds_number ** 0.9))) ** 2
    lerp = (reynolds_number - bot) / (top - bot)
    return lam * (1 - lerp) + turb * lerp

def relativeRoughness(roughness, diameter):
    return roughness / diameter

def estimateFrictionFactor(height, pipe_length, sin_theta, basin_area, pipe_area, pipe_diameter,
                            k_entry, kinematic_viscosity, k_junction=0.0, f_init=0.02):
    """
    Iteratively solves for friction factor.
    k_junction=0.0 for single pipe, set to junction loss coefficient for T-junction.
    V2 is velocity at junction inlet (pre-split).
    Continuity and energy are both satisfied by folding K_junction into the denominator.
    """
    coeff = pipe_area / basin_area
    f = f_init

    for _ in range(50):
        denominator = 1 - coeff**2 + f * (pipe_length / pipe_diameter) + k_entry + k_junction
        numerator = 2 * gravity * effectiveHead(height, pipe_diameter, pipe_length, sin_theta)

        if denominator <= 0 or numerator <= 0:
            break

        velocity = np.sqrt(numerator / denominator)
        reynolds = reynoldsNumber(velocity, pipe_diameter, kinematic_viscosity)
        rr = relativeRoughness(PVC_roughness, pipe_diameter)
        new_f = frictionFactor(rr, reynolds)

        if abs(new_f - f) < 1e-8:
            break
        f = new_f

    return f

def pipeChange(height, pipe_length, sin_theta, basin_area, pipe_area, pipe_diameter,
               k_entry, friction_factor, k_junction=0.0):
    """
    Computes dh/dt.
    For T-junction: V3 = V2/2 by continuity, but 2 branches means 2*A_p*V3 = A_p*V2.
    The factor of 2 from two branches cancels the /2 from continuity, so dh/dt = -(A_p/A_t)*V2.
    K_junction is folded into denominator to reduce V2, accounting for junction energy loss.
    """
    coeff = pipe_area / basin_area
    denominator = 1 - coeff**2 + friction_factor * (pipe_length / pipe_diameter) + k_entry + k_junction
    numerator = 2 * gravity * effectiveHead(height, pipe_diameter, pipe_length, sin_theta)

    if denominator <= 0 or numerator <= 0:
        return 0.0

    V2 = np.sqrt(numerator / denominator)
    return -coeff * V2

def rk4(height, pipe_length, sin_theta, basin_area, pipe_area, pipe_diameter,
        k_entry, kinematic_viscosity, k_junction=0.0, f_init=0.02, n=1000, s=0.01):

    def dhdt(h):
        f = estimateFrictionFactor(h, pipe_length, sin_theta, basin_area, pipe_area,
                                    pipe_diameter, k_entry, kinematic_viscosity, k_junction, f_init)
        return pipeChange(h, pipe_length, sin_theta, basin_area, pipe_area,
                          pipe_diameter, k_entry, f, k_junction)

    heights = [height]
    times = [0]
    curr_height = height
    curr_time = 0

    for _ in range(n):
        k1 = s * dhdt(curr_height)
        k2 = s * dhdt(curr_height + k1 / 2)
        k3 = s * dhdt(curr_height + k2 / 2)
        k4 = s * dhdt(curr_height + k3)

        curr_height += (k1 + 2*k2 + 2*k3 + k4) / 6
        curr_time += s

        heights.append(curr_height)
        times.append(curr_time)

        if height - curr_height > Config.drop_target:
            break

    return times, heights


def run_case(pipe_length, k_junction=0.0, label="Single pipe"):
    pipe_area = circleArea(Config.pipe_diameter)

    times, heights = rk4(
        height=Config.initial_height,
        pipe_length=pipe_length,
        sin_theta=Config.sin_theta,
        basin_area=Config.basin_area,
        pipe_area=pipe_area,
        pipe_diameter=Config.pipe_diameter,
        k_entry=Config.k_entry,
        kinematic_viscosity=kinematic_viscosity,
        k_junction=k_junction,
        n=Config.n_steps,
        s=Config.step_size
    )

    drop = Config.initial_height - heights[-1]
    reached = drop >= Config.drop_target
    return times[-1], heights[-1], reached


def main():
    pipe_lengths = [0.2, 0.3, 0.4, 0.6]

    print(f"{'L (cm)':<10} {'Single (s)':<14} {'T-junction (s)':<16} {'Diff (s)':<12} {'Diff (%)':<10}")
    print("-" * 62)

    for pipe_length in pipe_lengths:
        t_single, _, _ = run_case(pipe_length, k_junction=0.0)
        t_T, _, _ = run_case(pipe_length, k_junction=Config.K_junction)

        diff = t_T - t_single
        pct = diff / t_single * 100
        print(f"{pipe_length*100:<10.0f} {t_single:<14.1f} {t_T:<16.1f} {diff:<+12.1f} {pct:<+10.1f}")


if __name__ == "__main__":
    main()