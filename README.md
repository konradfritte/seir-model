# Extended SEIR-HD Model (with Reinfection)

A NumPy-based epidemiological simulation using the **Euler method**. This model extends the classic SEIR structure to include hospitalization and mortality dynamics, as well as the possibility of reinfection.

## 🚀 Features

- **Extended Compartments:** Tracks **S**usceptible, **E**xposed, **I**nfectious, **R**ecovered, **H**ospitalized, and **D**ead.
- **Reinfection Logic:** Models waning immunity where Recovered individuals can become Susceptible again.
- **Disease Profiles:** Includes pre-configured parameters for **COVID-19**, **Measles**, and **Influenza**.
- **NumPy Implementation:** Uses pure NumPy for iterative calculation via the Euler method.

## 🛠 Installation & Usage

Install the required dependencies:

```bash
pip install numpy matplotlib

Then run the simulation:
```bash
python main.py
