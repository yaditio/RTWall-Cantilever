# 🧱 RT Wall - Cantilever

**RT Wall - Cantilever** is a unified, interactive web application designed for the analysis, design, and verification of cantilever retaining walls. The application integrates classical geotechnical analytical calculations, soil stress distribution and bearing capacity analysis using **groundhog**, 2D Finite Element Method (FEM) simulation using **OpenSees**, ultimate concrete bending capacity checks using **concreteproperties**, and if necessary, soil-structure interaction via Winkler pile analysis using **openpile**.

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
- **Global Stability Analysis (2D SSRM)**: Evaluates the slope factor of safety (FS) using Finite Element Strength Reduction Method (SSRM) in `xslope` via the Griffiths & Lane (1999) displacement catastrophe sweep, displaying viscoplastic shear strain contour plots.
- **Pile Foundation Design**:
  - Distributes total retaining wall axial demand and overturning moments into a 2-row pile system.
  - Generates pile cross-section plots and verifies bending-axial interaction envelopes using `concrete-properties`.
  - Performs Winkler lateral response (deflection, shear, moment profiles) under load demands using `openpile`.
- **Interactive SVG Visualizations**: Proportionally scaled high-contrast drawing of the wall, soil layering, water tables, centroids, pile rows, and reaction forces.

---

## 🛠️ Local Installation & Run

### Prerequisites
- Python 3.8 to 3.10 is required (OpenPile and other dependencies do not support Python 3.11+).
- Cairo graphics library (needed by `drawsvg` for generating PNG/raster formats).

### Steps using Anaconda (Recommended)
If you do not have Anaconda installed, download it from the [Official Anaconda Website](https://www.anaconda.com/).

1. **Clone the repository**:
   ```bash
   git clone <your-repository-url>
   cd "RT Wall - Cantilever"
   ```

2. **Create and activate a Python 3.10 Conda environment**:
   ```bash
   conda create -n python310 python=3.10 -y
   conda activate python310
   ```

3. **Install the Cairo graphics dependency**:
   ```bash
   conda install -c conda-forge cairo -y
   ```
   *(Alternatively, on Ubuntu/Debian systems you can install it using: `sudo apt install libcairo2`)*

4. **Install Python dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

5. **Run the Streamlit app**:
   ```bash
   streamlit run App.py
   ```

### 📦 Standalone Desktop Executable (.exe)

The application can be compiled and run as a standalone, portable Windows executable (`.exe`), which bundles the Python environment, core libraries (including the OpenSees finite element solver), and Cairo rendering engine.

#### 1. Running the Executable
- Download the precompiled `RTWallCantilever.exe` from the **Releases** section of this repository.
- Double-click the downloaded executable to run it. A console window will open, start the local server, and automatically open your default web browser to the application page at `http://localhost:8501`.

#### 2. Re-compiling the Executable
If you modify `App.py` and wish to rebuild the portable executable:
1. Activate the Python 3.10 environment.
2. Install PyInstaller (if not already installed):
   ```bash
   pip install pyinstaller
   ```
3. Compile using the provided specification file:
   ```bash
   pyinstaller --clean --noconfirm RTWallCantilever.spec
   ```
   The output executable will be created/updated in the `dist/` folder.

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

## 📦 Core Libraries & Dependencies

This application uses the following libraries:
- [Streamlit](https://streamlit.io/) — Python web UI framework.
- [OpenSeesPy](https://openseespy.github.io/manual/index.html) — 2D/3D finite element analysis solver.
- [concreteproperties](https://concreteproperties.readthedocs.io/) — Concrete cross-section structural design verification.
- [sectionproperties](https://sectionproperties.readthedocs.io/) — Structural section property calculations.
- [openpile](https://github.com/pypile/openpile) — Winkler model implementation for lateral pile deflections.
- [xslope](https://xslope.org/) — Geotechnical slope stability and finite element seepage solver.
- [gmsh](https://gmsh.info/) — Three-dimensional finite element mesh generator.
- [groundhog](https://groundhog.readthedocs.io/) — Geotechnical engineering calculator library.
- [drawsvg](https://github.com/cduck/drawsvg) — SVG programmatic drawing generator.
- [NumPy](https://numpy.org/) & [pandas](https://pandas.pydata.org/) — Data parsing and numerical matrices.
- [Matplotlib](https://matplotlib.org/) — Geometry plotting and visualization engine.

---

## 📄 License

This project is licensed under the **GNU General Public License v3.0 (GPLv3)**. See the `LICENSE` file for details.
