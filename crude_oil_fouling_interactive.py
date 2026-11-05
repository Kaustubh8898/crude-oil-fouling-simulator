"""
Interactive crude-oil fouling model – student research prototype.
Run: python crude_oil_fouling_interactive.py
Optional dashboard: pip install -r requirements.txt && streamlit run app.py

Based on Coletti & Macchietto (2011) reference equations and geometry.
This is a reduced-order educational/research model, not plant-control software.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

RG = 8.314462618
DAY = 86400.0
L, DI, DO, NT, NP, NS = 6.1, 0.01905, 0.02540, 880, 4, 2
K_FRESH, K_COKED = 0.2, 1.0
EF, GAMMA, A_FOUL = 28491.0, 9.28e-13, 0.00165
RATE_SCALE_DEFAULT = 0.24906518298
EA = 50000.0
AGEING_RATES = {"slow": 1e-4, "intermediate": 1e-3, "fast": 1e-2}


def viscosity_cst_api(Tc, nu100):
    """API 11A4.1 estimate; guarded against nonphysical/out-of-range values."""
    Tf = Tc * 9 / 5 + 32
    B = np.log10(np.log10(nu100 + 0.7)) + 0.86960
    S = 0.28008 * np.log10(nu100) + 1.8616
    x = B * (559.67 / (Tf + 459.67)) ** S - 0.86960
    inner = 10 ** x
    nu = 10 ** inner - 0.7
    return float(np.clip(nu, 0.5, 1e5))


def fluid_props(Tc, fluid="crude"):
    """Approximate fluid properties used by this prototype."""
    nu100, api, cp = (13.9, 37.5, 2300.0) if fluid == "crude" else (28.8, 17.4, 2400.0)
    sg = 141.5 / (api + 131.5)
    rho_ref = sg * 999.016
    rho = max(rho_ref * (1 - 7e-4 * (Tc - 15.56)), 500.0)
    nu = viscosity_cst_api(Tc, nu100) * 1e-6
    mu = rho * nu
    Tk = Tc + 273.15
    k = max(0.08, 0.164 - 0.0001277 * Tk)
    return rho, mu, k, cp


def simulate(
    input_csv,
    duration_days=347,
    nodes=13,
    ageing="intermediate",
    rate_scale=RATE_SCALE_DEFAULT,
    shell_h=250.0,
    tube_flow_multiplier=1.0,
    shell_flow_multiplier=1.0,
):
    """Daily quasi-steady thermal solution + explicit daily deposit-state update."""
    data = pd.read_csv(input_csv)
    days = np.arange(0, min(duration_days, float(data.day.max())) + 1, 1, dtype=float)
    z = np.linspace(0, L, nodes)
    dz = L / (nodes - 1)
    weights = np.ones(nodes)
    weights[[0, -1]] = 0.5
    Rf = np.zeros((NS, NP, nodes))
    delta = np.full_like(Rf, 1e-7)
    youth = np.ones_like(Rf)
    hist, local = [], []

    for day in days:
        f = {c: float(np.interp(day, data.day, data[c])) for c in data.columns if c != "day"}
        Qt = f["tube_flow_m3_h"] / 3600 * tube_flow_multiplier
        Qs = f["shell_flow_m3_h"] / 3600 * shell_flow_multiplier
        Tin = f["crude_inlet_C"]
        Ts_in = f["residuum_inlet_C"]
        Ts_out = f["residuum_outlet_C"]
        rho_t, mu_t, k_t, cp_t = fluid_props(Tin, "crude")
        rho_s, mu_s, k_s, cp_s = fluid_props((Ts_in + Ts_out) / 2, "residuum")
        mdot_t = rho_t * Qt
        mdot_s = rho_s * Qs
        tube_bulk = np.full((NS, NP, nodes), Tin)
        tube_bulk[0, :, 0] = Tin
        q_local = np.zeros_like(tube_bulk)
        hi_local = np.zeros_like(tube_bulk)
        re_local = np.zeros_like(tube_bulk)
        tau_local = np.zeros_like(tube_bulk)
        shell_bulk_in = Ts_in
        total_dp = 0.0
        total_tube_duty = 0.0

        for shell in range(NS):
            shell_direction = 1 if shell == 0 else -1
            shell_bulk = np.linspace(shell_bulk_in, Ts_out if shell == 1 else (Ts_in + Ts_out) / 2, nodes)
            for p in range(NP):
                for j in range(nodes):
                    Tb = tube_bulk[shell, p, j]
                    rho, mu, k, cp = fluid_props(Tb, "crude")
                    rflow = max(DI / 2 - delta[shell, p, j], DI * 0.15)
                    area = np.pi * rflow**2 * NT / NP
                    velocity = Qt / max(area, 1e-12)
                    Re = rho * velocity * (2 * rflow) / max(mu, 1e-12)
                    Pr = cp * mu / max(k, 1e-12)
                    Nu = 0.023 * max(Re, 1.0) ** 0.8 * max(Pr, 0.1) ** 0.4
                    hi = k * Nu / max(2 * rflow, 1e-12)
                    cf = 0.0035 + 0.264 * max(Re, 1.0) ** -0.42
                    tau = cf * rho * velocity**2 / 2
                    hs = shell_h
                    U = 1 / (1 / max(hi, 1e-9) + Rf[shell, p, j] + 1 / max(hs, 1e-9))
                    q = U * np.pi * DO * dz * (shell_bulk[j] - Tb)
                    q_local[shell, p, j] = q
                    hi_local[shell, p, j] = hi
                    re_local[shell, p, j] = Re
                    tau_local[shell, p, j] = tau
                    total_dp += cf * rho * velocity**2 / max(rflow, 1e-9) * dz
                    total_tube_duty += q * NT / NP

                # Approximate passwise tube heating; flow reverses direction between passes.
                if p < NP - 1:
                    pass_duty = np.sum(q_local[shell, p, :] * weights) * NT / NP * dz
                    dT = pass_duty / max(mdot_t * cp_t, 1e-9)
                    tube_bulk[shell, p + 1, :] = np.linspace(Tb + dT, Tb + dT, nodes)

            # Advance shell-side bulk temperature using the local exchanged duty.
            shell_duty = np.sum(q_local[shell] * weights[None, :]) * NT / NP * dz
            shell_bulk_in -= shell_direction * shell_duty / max(mdot_s * cp_s, 1e-9)

        # Deposit kinetics and ageing update.
        for shell in range(NS):
            for p in range(NP):
                for j in range(nodes):
                    Tb = tube_bulk[shell, p, j]
                    Tf = Tb + 0.55 * max(shell_bulk_in - Tb, 0.0)
                    rho, mu, k, cp = fluid_props(Tb, "crude")
                    rflow = max(DI / 2 - delta[shell, p, j], DI * 0.15)
                    area = np.pi * rflow**2 * NT / NP
                    velocity = Qt / max(area, 1e-12)
                    Re = rho * velocity * (2 * rflow) / max(mu, 1e-12)
                    Pr = cp * mu / max(k, 1e-12)
                    tau = tau_local[shell, p, j]
                    formation = rate_scale * A_FOUL * max(Re, 1.0) ** -0.66 * max(Pr, 0.1) ** -0.33
                    formation *= np.exp(-EF / (RG * max(Tf + 273.15, 1.0)))
                    erosion = GAMMA * tau
                    dRf = max(formation - erosion, 0.0) * DAY
                    Rf[shell, p, j] += dRf
                    delta[shell, p, j] += K_FRESH * dRf
                    youth[shell, p, j] *= np.exp(-AGEING_RATES[ageing] * DAY)
                    youth[shell, p, j] = np.clip(youth[shell, p, j], 0.0, 1.0)

                    local.append({
                        "day": day, "shell": shell + 1, "pass": p + 1, "z_m": z[j],
                        "deposit_m": delta[shell, p, j], "Rf_m2K_W": Rf[shell, p, j],
                        "film_C": Tf, "Re": Re, "Pr": Pr, "wall_shear_Pa": tau,
                        "deposit_youth": youth[shell, p, j],
                    })

        mean_delta = float(np.mean(delta))
        max_delta = float(np.max(delta))
        tube_duty_MW = total_tube_duty / 1e6
        shell_duty_MW = mdot_s * cp_s * (Ts_in - Ts_out) / 1e6
        tube_outlet = Tin + tube_duty_MW * 1e6 / max(mdot_t * cp_t, 1e-9)
        hist.append({
            "day": day,
            "mean_deposit_thickness_mm": mean_delta * 1e3,
            "max_deposit_thickness_mm": max_delta * 1e3,
            "raw_pressure_drop_bar": total_dp / 1e5,
            "tube_outlet_model_C": tube_outlet,
            "tube_duty_MW": tube_duty_MW,
            "shell_duty_from_inputs_MW": shell_duty_MW,
        })

    return pd.DataFrame(hist), pd.DataFrame(local)


def main():
    parser = argparse.ArgumentParser(description="Interactive crude-oil fouling model prototype")
    parser.add_argument("--days", type=int, default=347)
    parser.add_argument("--nodes", type=int, default=13, choices=[7, 13, 25, 49])
    parser.add_argument("--ageing", choices=list(AGEING_RATES), default="intermediate")
    parser.add_argument("--input", default=str(Path(__file__).with_name("Figure8_digitised_daily_inputs.csv")))
    args = parser.parse_args()
    out_dir = Path("model_output")
    out_dir.mkdir(exist_ok=True)
    summary, local = simulate(args.input, args.days, args.nodes, args.ageing)
    summary.to_csv(out_dir / "daily_summary.csv", index=False)
    local.to_csv(out_dir / "local_axial_results.csv", index=False)
    print(summary.tail(1).to_string(index=False))
    print(f"\nWrote results to: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
