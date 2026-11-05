# Crude-oil fouling interactive simulator

A runnable, reduced-order research prototype based on the reference case in Coletti & Macchietto (2011). The Streamlit dashboard recalculates model outputs when scenario inputs are changed.

> **Research prototype only.** This is not connected to live refinery sensors and must not be used for plant operation, equipment protection, or safety decisions. Energy closure and physical validation are not established.

## Contents

- `crude_oil_fouling_interactive.py` — simulation engine and command-line runner.
- `app.py` — interactive Streamlit dashboard.
- `Figure8_digitised_daily_inputs.csv` — approximate daily inputs digitised visually from Figure 8 of the reference paper.
- `requirements.txt` — Python dependencies.

## Run the dashboard locally

Requires Python 3.10 or newer.

```bash
python -m venv .venv
# macOS / Linux
source .venv/bin/activate
# Windows PowerShell
# .venv\Scripts\Activate.ps1

pip install -r requirements.txt
streamlit run app.py
```

Open the local URL printed by Streamlit. Change the controls in the sidebar and select **Run / update simulation**. The dashboard displays time trends, axial profiles, model estimates, and CSV downloads.

## Run from the command line

```bash
python crude_oil_fouling_interactive.py --days 347 --nodes 13 --ageing intermediate
```

Outputs are written to `model_output/`.

## Model controls

- Simulation duration: 30–347 days.
- Axial nodes per pass: 7, 13, 25, or 49.
- Deposit ageing scenario: slow, intermediate, or fast.
- Fouling-rate multiplier.
- Shell-side heat-transfer coefficient.
- Tube-side and shell-side flow multipliers.

## Scientific limitations

- Figure 8 input profiles are visual digitisation, not original refinery historian data.
- The shell-side heat-transfer coefficient is assumed.
- Heat capacities are fixed; density and thermal-conductivity relations are approximate.
- The fouling-rate multiplier is provisional.
- This is a simplified one-dimensional model. Results require energy-balance review, calibration, and physical validation against the reference figures before being treated as predictive engineering results.
- Scenario recalculation is interactive; it is not a live plant-data feed.

## Reference

Coletti, F. & Macchietto, S. (2011). “A Dynamic, Distributed Model of Shell-and-Tube Heat Exchangers Undergoing Crude Oil Fouling.” *Industrial & Engineering Chemistry Research*, 50, 4515–4533. https://doi.org/10.1021/ie901991g
