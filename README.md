# 🧱 RT Wall - Cantilever

**RT Wall - Cantilever** is a unified, interactive web application designed for the analysis, design, and verification of cantilever retaining walls and pile foundations. The application integrates classical geotechnical analytical calculations, 2D Finite Element Method (FEM) simulation using **OpenSees**, ultimate concrete bending capacity checks using **concreteproperties**, and soil-structure interaction via Winkler pile analysis using **openpile**.

---

## ⚠️ Disclaimer

This project is currently experimental and under active development.

It is provided for educational and learning purposes. Features may change without notice, and bugs or unexpected behavior may occur.

This software is not a certified engineering design tool and should not be used as the sole basis for engineering decisions. Always validate results independently and comply with applicable codes, standards, and professional engineering practices.

---

## 🚀 Features

- **Geometric Design**: Custom wall height, toe/heel slab sizing, tapering, and soil cover inputs.
- **Geotechnical Surcharge**: Integrates strip loads or uniform surcharges using Boussinesq distributions.
- **Seismic Analysis**: Automatically applies pseudo-static seismic loads ($k_h = \text{PGA} \times \text{FPGA}$) to OpenSees nodes.
- **Concrete Strength Verification**: Section capacity check for Stem Wall, Footing Toe, and Footing Heel.
- **Pile Foundation Design**:
  - Distributes total retaining wall axial demand and overturning moments into a 2-row pile system.
  - Generates pile cross-section plots and verifies bending-axial interaction envelopes using `concrete-properties`.
  - Performs Winkler lateral response (deflection, shear, moment profiles) under load demands using `openpile`.
- **Interactive SVG Visualizations**: Proportionally scaled high-contrast drawing of the wall, soil layering, water tables, centroids, pile rows, and reaction forces.

---

## 🛠️ Local Installation & Run

### Prerequisites
- Python 3.8 to 3.10 (OpenPile and dependencies do not currently support Python 3.11+).
- Cairo graphics library (needed by `drawsvg` for PNG/raster image outputs). Install on Ubuntu/Debian via:
  ```bash
  sudo apt install libcairo2
  ```

### Steps
1. Clone the repository and navigate into the folder:
   ```bash
   git clone <your-repository-url>
   cd "RT Wall - Cantilever"
   ```
2. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Run the Streamlit web server:
   ```bash
   streamlit run App.py
   ```

---

## 🐳 Docker Deployment

The application is prepared for deployment inside a Docker container.

### Build the Docker Image
```bash
docker build -t rt-wall-cantilever .
```

### Run the Container
```bash
docker run -d -p 8501:8501 rt-wall-cantilever
```
Once started, access the app at `http://localhost:8501`.

---

## 📄 License

This project is licensed under the **GNU General Public License v3.0 (GPLv3)**. See the `LICENSE` file for details.
