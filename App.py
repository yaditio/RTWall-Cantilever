# ----------------------------------------
# Cantilever Retaining Wall
# ----------------------------------------

# Importing necessary libraries
import openseespy.opensees as ops
import opsvis as opsv
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
import drawsvg as draw
import pandas as pd

# Local implementations of Boussinesq strip load stress distribution functions (ensures Groundhog version compatibility)
def stresses_stripload_local(z, x, width, imposedstress, triangular=False):
    z = max(1e-5, float(z))
    R_1 = np.sqrt(x ** 2 + z ** 2)
    R_2 = np.sqrt((x - width) ** 2 + z ** 2)

    _theta1 = np.arccos(z / R_1)
    _theta2 = np.arccos(z / R_2)
    if x < width:
        _theta2 = -_theta2
    beta = _theta2
    alpha = _theta1 - beta

    if triangular:
        _delta_sigma_z = (imposedstress / np.pi) * (
            (x / width) * alpha -
            0.5 * np.sin(2 * beta)
        )
        _delta_sigma_x = (imposedstress / np.pi) * (
            (x / width) * alpha -
            (z / width) * np.log((R_1 ** 2) / (R_2 ** 2)) +
            0.5 * np.sin(2 * beta)
        )
        _delta_tau_zx = (imposedstress / (2 * np.pi)) * (
            1 +
            np.cos(2 * beta) -
            2 * (z / width) * alpha
        )
    else:
        _delta_sigma_z = (imposedstress / np.pi) * (
            alpha +
            np.sin(alpha) * np.cos(alpha + 2 * beta)
        )
        _delta_sigma_x = (imposedstress / np.pi) * (
            alpha -
            np.sin(alpha) * np.cos(alpha + 2 * beta)
        )
        _delta_tau_zx = (imposedstress / np.pi) * (
            np.sin(alpha) * np.sin(alpha + 2 * beta)
        )

    return {
        'delta sigma z [kPa]': _delta_sigma_z,
        'delta sigma x [kPa]': _delta_sigma_x,
        'delta tau zx [kPa]': _delta_tau_zx,
    }

def stresses_stripload_retainingwall_local(imposedstress, width, offset, toe_depth, depth):
    H0 = max(1e-5, float(toe_depth))
    z = max(1e-5, float(depth))
    a = float(offset)
    B = float(width)
    
    alpha_rad = np.arctan(a / z)
    alpha_far_rad = np.arctan((a + B) / z)
    beta_rad = alpha_far_rad - alpha_rad
    alpha_bisector_rad = alpha_rad + beta_rad / 2.0
    
    _delta_sigma_x = (2.0 * imposedstress / np.pi) * (
        beta_rad - np.sin(beta_rad) * np.cos(2.0 * alpha_bisector_rad)
    )
    
    theta_1 = np.degrees(np.arctan(a / H0))
    theta_2 = np.degrees(np.arctan((a + B) / H0))
    R_1 = ((a + B) ** 2) * (90.0 - theta_2)
    R_2 = (a ** 2) * (90.0 - theta_1)
    
    _delta_P_x = (imposedstress / 90.0) * (H0 * (theta_2 - theta_1))
    
    denominator = 2.0 * H0 * (theta_2 - theta_1)
    if abs(denominator) < 1e-5:
        _z_bar = H0 / 2.0
    else:
        _z_bar = (H0**2 * (theta_2 - theta_1) - (R_1 - R_2) + 57.3 * B * H0) / denominator
        
    return {
        'delta sigma x [kPa]': _delta_sigma_x,
        'delta P x [kN/m]': _delta_P_x,
        'z bar [m]': _z_bar,
        'theta 1 [deg]': theta_1,
        'theta 2 [deg]': theta_2,
        'R 1 [m]': R_1,
        'R 2 [m]': R_2,
        'alpha [deg]': np.degrees(alpha_bisector_rad),
        'beta [deg]': np.degrees(beta_rad)
    }

# Page Configuration for wide layout
st.set_page_config(layout="wide", page_title="Cantilever Retaining Wall Studio", page_icon="🧱")

# CSS Injection for Premium Styling
st.markdown("""
<style>
    /* Gradient Dashboard Header Banner */
    .dashboard-header {
        background: linear-gradient(135deg, #1f4068 0%, #162447 100%);
        color: #ffffff;
        padding: 24px;
        border-radius: 12px;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.15);
        margin-bottom: 24px;
        text-align: center;
    }
    .dashboard-header h1 {
        margin: 0;
        font-family: 'Outfit', 'Inter', sans-serif;
        font-size: 2.5rem;
        font-weight: 700;
        letter-spacing: -0.5px;
    }
    .dashboard-header p {
        margin: 8px 0 0 0;
        font-size: 1.1rem;
        opacity: 0.9;
        font-weight: 300;
    }
    /* Section subtitle details */
    .section-desc {
        color: #6c757d;
        font-size: 0.9rem;
        margin-bottom: 12px;
    }
</style>
""", unsafe_allow_html=True)

# ----------------------------------------

# Clean up anything created by invoked opensees commands
ops.wipe()

# We define the space in which we will create the model
ops.model('basic','-ndm',2,'-ndf',2) # Displacements in 2 directions, out-of-plane displacements and rotations are restricted.

# Streamlit sidebar inputs (re-runs app on change)
tab_geom, tab_soil, tab_struct, tab_seismic, tab_reinf, tab_bearing, tab_pile = st.sidebar.tabs(['Geometry', 'Soil', 'Structure', 'Seismic', 'Concrete & Reinf', 'Bearing Capacity', 'Pile'])

with tab_geom:
    t = st.number_input('Out-of-plane thickness (m)', value=1.0, step=0.01, format="%.3f")
    Hw = st.number_input('Wall height (m)', value=4.5, step=0.01, format="%.3f")
    toe = st.number_input('Toe length (m)', value=1.0, step=0.01, format="%.3f")
    heel = st.number_input('Heel length (m)', value=1.75, step=0.01, format="%.3f")
    top_wall = st.number_input('Top wall thickness (m)', value=0.25, step=0.01, format="%.3f")
    bot_wall = st.number_input('Bottom wall thickness (m)', value=0.45, step=0.01, format="%.3f")
    h_ftg = st.number_input('Footing thickness (m)', value=0.25, step=0.01, format="%.3f")

with tab_soil:
    # Soil parameters
    gamma_soil_dry = st.number_input('Dry soil unit weight (kN/m3)', value=18.0, step=0.1, format="%.2f")
    gamma_soil_wet = st.number_input('Wet soil unit weight (kN/m3)', value=20.0, step=0.1, format="%.2f")
    phi = st.number_input('Friction angle (deg)', value=45, step=1)
    surcharge_type = st.selectbox('Surcharge Type', ['Uniform', 'Strip Load'], index=0)
    q = st.number_input('Surcharge load q (kN/m2)', value=10.0, step=0.1, format="%.2f")
    if surcharge_type == 'Strip Load':
        width_surcharge = st.number_input('Surcharge width B (m)', value=2.0, step=0.1, format="%.2f")
        offset_surcharge = st.number_input('Surcharge offset from wall face a (m)', value=1.5, step=0.1, format="%.2f")
    else:
        width_surcharge = 1e5
        offset_surcharge = 0.0
    h_soil = st.number_input('Soil height above heel (m)', value=4.5, step=0.01, format="%.3f")
    h_soil_toe = st.number_input('Soil above toe (m)', value=0.1, step=0.01, format="%.3f")

    # Water parameters
    gamma_w = st.number_input('Water unit weight (kN/m3)', value=9.81, step=0.01, format="%.2f")
    Hwtr = st.number_input('Water height behind wall (m)', value=0.5, step=0.01, format="%.2f")
    Hwtr_front = st.number_input('Water height in front of wall (m)', value=0.15, step=0.01, format="%.3f")

with tab_seismic:
    #Seismic Parameter
    PGA = st.number_input('Peak Ground Acceleration (g)', value=0.4123, step=0.01, format="%.2f")
    FPGA = st.number_input('Seismic Load Factor', value=1.2, step=0.01, format="%.2f")

with tab_struct:
    # Concrete self-weight
    gamma_c = st.number_input('Unit weight of concrete (kN/m3)', value=24.0, step=0.1, format="%.2f")

with tab_reinf:
    st.subheader("Concrete Properties")
    fc = st.number_input('Comp. strength f\'c (MPa)', value=50.0, step=5.0)
    Ec = st.number_input('Elastic Modulus Ec (MPa)', value=34.8e3, step=1000.0)
    ft_tensile = st.number_input('Flexural tensile strength (MPa)', value=4.2, step=0.1)
    
    st.subheader("Steel Properties")
    fy = st.number_input('Yield strength fsy (MPa)', value=500.0, step=50.0)
    Es = st.number_input('Elastic Modulus Es (MPa)', value=200e3, step=5000.0)
    
    # Indonesian Standard SNI Rebar diameter lists in mm
    sni_diameters = [6, 8, 10, 12, 13, 14, 16, 19, 22, 25, 29, 32, 36, 40, 50]
    
    st.subheader("Stem Reinforcement")
    stem_rebar_dia = st.selectbox('Stem main bar diameter (mm)', sni_diameters, index=9, help="SNI diameter options (e.g. 25 = D25)")
    stem_bar_area = np.pi * (stem_rebar_dia ** 2) / 4.0
    st.caption(f"Calculated Stem Main Bar Area: **{stem_bar_area:.2f} mm²**")
    
    stem_spacing_x = st.number_input('Stem bar spacing x (mm)', value=150.0, step=10.0)
    
    stem_sidebar_dia = st.selectbox('Stem side bar diameter (mm)', sni_diameters, index=6)
    stem_sidebar_area = np.pi * (stem_sidebar_dia ** 2) / 4.0
    st.caption(f"Calculated Stem Side Bar Area: **{stem_sidebar_area:.2f} mm²**")
    
    stem_spacing_y = st.number_input('Stem side bar spacing y (mm)', value=100.0, step=10.0)
    stem_cover = st.number_input('Stem concrete cover (mm)', value=50.0, step=5.0)
    
    st.subheader("Footing Reinforcement")
    ftg_rebar_dia = st.selectbox('Footing main bar diameter (mm)', sni_diameters, index=9)
    ftg_bar_area = np.pi * (ftg_rebar_dia ** 2) / 4.0
    st.caption(f"Calculated Footing Main Bar Area: **{ftg_bar_area:.2f} mm²**")
    
    ftg_spacing_x = st.number_input('Footing bar spacing x (mm)', value=150.0, step=10.0)
    
    ftg_sidebar_dia = st.selectbox('Footing side bar diameter (mm)', sni_diameters, index=6)
    ftg_sidebar_area = np.pi * (ftg_sidebar_dia ** 2) / 4.0
    st.caption(f"Calculated Footing Side Bar Area: **{ftg_sidebar_area:.2f} mm²**")
    
    ftg_spacing_y = st.number_input('Footing side bar spacing y (mm)', value=200.0, step=10.0)
    ftg_cover = st.number_input('Footing concrete cover (mm)', value=50.0, step=5.0)

with tab_bearing:
    enable_bearing = st.checkbox('Enable Bearing Capacity Analysis', value=False)
    if enable_bearing:
        bearing_soil_type = st.selectbox('Soil Type for Bearing', ['Sand (Drained)', 'Clay (Undrained)'], index=0)
        if bearing_soil_type == 'Clay (Undrained)':
            su_val = st.number_input('Undrained shear strength su (kPa)', value=50.0, step=5.0, format="%.1f")
        else:
            su_val = 0.0
        FS_bearing = st.number_input('Bearing Capacity Safety Factor (FS)', value=3.0, step=0.1, format="%.2f")
    else:
        bearing_soil_type = 'Sand (Drained)'
        su_val = 0.0
        
with tab_pile:
    enable_pile = st.checkbox('Enable Pile Foundation Analysis', value=False)
    if enable_pile:
        pile_shape = st.selectbox('Pile Shape', ['Circle', 'Rectangle'], index=0)
        pile_material = st.selectbox('Pile Material', ['Concrete', 'Steel'], index=0)
        pile_offset = st.number_input('Offset from Edge to Pile Centroid (m)', value=0.50, step=0.05, format="%.2f")
        L_pile = st.number_input('Pile Length (m)', value=10.0, step=0.5, format="%.1f")
        
        if pile_shape == 'Circle':
            diameter_pile = st.number_input('Pile Diameter (m)', value=0.40, step=0.05, format="%.2f")
            width_x_pile = diameter_pile
            width_y_pile = diameter_pile
        else:
            width_x_pile = st.number_input('Pile Width in Bending direction dx (m)', value=0.40, step=0.05, format="%.2f")
            width_y_pile = st.number_input('Pile Width perpendicular to Bending dy (m)', value=0.40, step=0.05, format="%.2f")
            diameter_pile = width_x_pile
            
        if pile_material == 'Concrete':
            n_rebar_pile = st.number_input('Number of Pile Rebars', value=8, min_value=4, step=1)
            rebar_dia_pile = st.selectbox('Pile Rebar Diameter (mm)', sni_diameters, index=9)
            cover_pile = st.number_input('Pile Concrete Cover (mm)', value=50.0, step=5.0)
        else:
            n_rebar_pile = 8
            rebar_dia_pile = 25
            cover_pile = 50.0
        
        st.subheader("Soil Profile under Footing")
        pile_soil_type = st.selectbox('Soil Type for Piles', ['Sand (Drained)', 'Clay (Undrained)'], index=0)
        if pile_soil_type == 'Sand (Drained)':
            phi_pile_soil = st.number_input('Pile Soil Friction Angle (deg)', value=30.0, step=1.0, format="%.1f")
            cohesion_pile_soil = st.number_input('Pile Soil Cohesion (kPa)', value=0.0, step=1.0, format="%.1f")
            gamma_pile_soil = st.number_input('Pile Soil Unit Weight (kN/m3)', value=18.0, step=0.5, format="%.1f")
            su_pile_soil = 0.0
            eps50_pile_soil = 0.015
        else:
            su_pile_soil = st.number_input('Pile Soil Undrained Shear Strength su (kPa)', value=50.0, step=5.0, format="%.1f")
            gamma_pile_soil = st.number_input('Pile Soil Unit Weight (kN/m3)', value=18.0, step=0.5, format="%.1f")
            eps50_pile_soil = st.number_input('Pile Soil Strain at 50% max stress eps50', value=0.015, step=0.005, format="%.3f")
            phi_pile_soil = 30.0
            cohesion_pile_soil = 0.0
            
        st.subheader("Analysis Safety Factors")
        FS_pile_axial = st.number_input('Axial Safety Factor (piles)', value=2.5, step=0.1, format="%.2f")
        FS_pile_lateral = st.number_input('Lateral Safety Factor (piles)', value=2.5, step=0.1, format="%.2f")
        pile_loading_type = st.selectbox('OpenPile Loading Type', ['static', 'cyclic'], index=0)
    else:
        pile_shape = 'Circle'
        pile_material = 'Concrete'
        pile_offset = 0.50
        L_pile = 10.0
        diameter_pile = 0.40
        width_x_pile = 0.40
        width_y_pile = 0.40
        n_rebar_pile = 8
        rebar_dia_pile = 25
        cover_pile = 50.0
        pile_soil_type = 'Sand (Drained)'
        phi_pile_soil = 30.0
        cohesion_pile_soil = 0.0
        su_pile_soil = 50.0
        gamma_pile_soil = 18.0
        eps50_pile_soil = 0.015
        FS_pile_axial = 2.5
        FS_pile_lateral = 2.5
        pile_loading_type = 'static'

# Derived geometry
ftg = toe + heel + top_wall + bot_wall  # total footing length
taper = bot_wall - top_wall
taper_length = np.sqrt(taper**2 + Hw**2)

# Build a refined mesh by splitting the original triangles while keeping
# the overall outer geometry unchanged. We create 8 nodes along the
# footing (bottom and top of footing) and 8 stem layers (back/front
# pairs) resulting in 32 nodes and 30 Tri31 elements.

# compute x positions along footing by subdividing the original three
# segments (0->toe, toe->toe+bot_wall, toe+bot_wall->ftg) into 3,2,2
# subdivisions respectively so total nodes along footing = 1+3+2+2 = 8
seg1_div = 3
seg2_div = 2
seg3_div = 2

xs = []
# Segment 1: 0 -> toe
for i in range(seg1_div):
    xs.append(0.0 + (toe - 0.0) * (i / float(seg1_div)))
xs.append(float(toe))
# Segment 2: toe -> toe+bot_wall
for i in range(1, seg2_div + 1):
    xs.append(float(toe) + (bot_wall) * (i / float(seg2_div)))
# Segment 3: toe+bot_wall -> ftg
for i in range(1, seg3_div + 1):
    xs.append((toe + bot_wall) + (ftg - (toe + bot_wall)) * (i / float(seg3_div)))

# ensure unique and sorted (floating tolerance)
xs_unique = []
for x in xs:
    if not any(abs(x - xx) < 1e-9 for xx in xs_unique):
        xs_unique.append(x)

# bottom nodes (y=0): IDs 1..8
num_foot_nodes = len(xs_unique)
bottom_node_ids = list(range(1, num_foot_nodes + 1))
for idx, x in enumerate(xs_unique):
    nid = bottom_node_ids[idx]
    ops.node(nid, float(x), 0.0)

# top-of-footing nodes (y=h_ftg): IDs 9..16
top_node_ids = list(range(num_foot_nodes + 1, num_foot_nodes * 2 + 1))
for idx, x in enumerate(xs_unique):
    nid = top_node_ids[idx]
    ops.node(nid, float(x), h_ftg)

# Stem nodes: create 8 back/front pairs above the top-of-footing
# back nodes follow the tapered profile, front nodes at x = toe+bot_wall
# IDs 17..32 (back1,front1, back2,front2, ...)
stem_pairs = []
for i in range(1, 9):
    z = h_ftg + Hw * (i / 8.0)
    back_x = toe + taper * (i / 8.0)
    front_x = toe + bot_wall
    back_id = num_foot_nodes * 2 + (2 * i - 1)
    front_id = num_foot_nodes * 2 + (2 * i)
    ops.node(back_id, float(back_x), float(z))
    ops.node(front_id, float(front_x), float(z))
    stem_pairs.append((back_id, front_id))

# Material definition
ops.nDMaterial('ElasticIsotropic', 1, 32000000.0, 0.2)

# Element definition
# Footing: create two triangles per cell between bottom/top nodes
elem_id = 1
for i in range(num_foot_nodes - 1):
    b1 = bottom_node_ids[i]
    b2 = bottom_node_ids[i + 1]
    t1 = top_node_ids[i]
    t2 = top_node_ids[i + 1]
    # triangle 1: b1, t2, t1
    ops.element('Tri31', elem_id, b1, t2, t1, 1.0, 'PlaneStrain', 1)
    elem_id += 1
    # triangle 2: b1, b2, t2
    ops.element('Tri31', elem_id, b1, b2, t2, 1.0, 'PlaneStrain', 1)
    elem_id += 1

# Stem: connect top-of-footing (left_col, right_col) to stem pairs
# find indices of left (x==toe) and front (x==toe+bot_wall) among top_node x positions
top_xs = xs_unique
def find_index_close(val, arr):
    for ii, xx in enumerate(arr):
        if abs(xx - val) < 1e-6:
            return ii
    return None

left_col_idx = find_index_close(toe, top_xs)
front_col_idx = find_index_close(toe + bot_wall, top_xs)
if left_col_idx is None or front_col_idx is None:
    # fallback to using the two central columns if exact matches fail
    left_col_idx = 1
    front_col_idx = 2

left_node = top_node_ids[left_col_idx]
right_node = top_node_ids[front_col_idx]

# build stem elements following the splitting pattern used previously
stem_elements_by_layer = []
for j in range(1, 9):
    b_j, f_j = stem_pairs[j - 1]
    if j == 1:
        # first strip connects left_node, front_j, back_j and right_node, front_j, left_node
        ops.element('Tri31', elem_id, left_node, f_j, b_j, 1.0, 'PlaneStrain', 1)
        e1 = elem_id
        elem_id += 1
        ops.element('Tri31', elem_id, right_node, f_j, left_node, 1.0, 'PlaneStrain', 1)
        e2 = elem_id
        elem_id += 1
    else:
        b_prev, f_prev = stem_pairs[j - 2]
        # element connecting previous back to current back via current front
        ops.element('Tri31', elem_id, b_prev, f_j, b_j, 1.0, 'PlaneStrain', 1)
        e1 = elem_id
        elem_id += 1
        # element connecting previous front, current front, previous back
        ops.element('Tri31', elem_id, f_prev, f_j, b_prev, 1.0, 'PlaneStrain', 1)
        e2 = elem_id
        elem_id += 1
    stem_elements_by_layer.append([e1, e2])

# total nodes and elements
total_nodes = num_foot_nodes * 2 + len(stem_pairs) * 2
total_elements = elem_id - 1

# heel_nodes kept for reference (top-of-footing sampling points)
heel_nodes = top_node_ids

st.success('Parameters updated — Streamlit will re-run automatically when inputs change.')

# Defining time series
ops.timeSeries("Linear", 1)

# Solution algorithm
ops.algorithm('Linear')

# Define a load class
# We have only one type of load.
ops.pattern("Plain", 1, 1)

# Active earth pressure coefficient (Rankine)
Ka = np.tan(np.radians(45 - phi/2))**2

# Back face wall nodes (use front nodes of stem pairs as sampling points)
wall_nodes = {}
for i, (b, f) in enumerate(stem_pairs, start=1):
    wall_nodes[f] = Hw * (i / 8.0)

def polygon_centroid(vertices):
    # vertices: list of (x,y) tuples (must be closed or will close)
    A = 0.0
    Cx = 0.0
    Cy = 0.0
    n = len(vertices)
    for i in range(n):
        x0, y0 = vertices[i]
        x1, y1 = vertices[(i + 1) % n]
        cross = x0 * y1 - x1 * y0
        A += cross
        Cx += (x0 + x1) * cross
        Cy += (y0 + y1) * cross
    A *= 0.5
    if abs(A) < 1e-12:
        return (0.0, 0.0)
    Cx /= (6.0 * A)
    Cy /= (6.0 * A)
    return (Cx, Cy)

# Build stem polygon (geometry of the wall above top-of-footing)
stem_poly = []
try:
    stem_poly.append((toe, -h_ftg))
    stem_poly.append((toe + taper/4.0, -(h_ftg + Hw/4.0)))
    stem_poly.append((toe + taper/2.0, -(h_ftg + Hw/2.0)))
    stem_poly.append((toe + (taper * 3.0 / 4.0), -(h_ftg + Hw * 3.0 / 4.0)))
    stem_poly.append((toe + taper, -(h_ftg + Hw)))
    stem_poly.append((toe + bot_wall, -(h_ftg + Hw)))
    stem_poly.append((toe + bot_wall, -h_ftg))
except Exception:
    stem_poly = []

stem_centroid = (None, None)
if stem_poly:
    stem_centroid = polygon_centroid(stem_poly)

# Base (footing slab) polygon (rectangular approximation)
base_poly = [(0.0, 0.0), (float(ftg), 0.0), (float(ftg), -h_ftg), (0.0, -h_ftg)]
base_centroid = polygon_centroid(base_poly)

# Soil (heel block) polygons for drawing and centroids
h_wet_height = min(h_soil, Hwtr)
soil_dry_poly = [(toe + bot_wall, -(h_ftg + h_wet_height)), (float(ftg), -(h_ftg + h_wet_height)), (float(ftg), -(h_ftg + h_soil)), (toe + bot_wall, -(h_ftg + h_soil))] if h_soil > Hwtr else []
soil_dry_centroid = polygon_centroid(soil_dry_poly) if soil_dry_poly else (None, None)

soil_wet_poly = [(toe + bot_wall, -h_ftg), (float(ftg), -h_ftg), (float(ftg), -(h_ftg + h_wet_height)), (toe + bot_wall, -(h_ftg + h_wet_height))] if h_wet_height > 0 else []
soil_wet_centroid = polygon_centroid(soil_wet_poly) if soil_wet_poly else (None, None)

# Map each stem front node to stem centroid x for previous uses
centroid_x_wall = {}
for _, front_id in stem_pairs:
    centroid_x_wall[front_id] = stem_centroid[0] if stem_centroid[0] is not None else (toe + bot_wall) / 2.0

# Compute lateral earth pressure at each node
sigma_total = {}
sigma_soil_q = {}
sigma_w_dict = {}

for node, z in wall_nodes.items():

    # Soil pressure (z is elevation from top of footing)
    h_dry_z = max(0.0, h_soil - max(z, Hwtr))
    h_wet_z = max(0.0, min(h_soil, Hwtr) - z)
    sigma_soil = Ka * (gamma_soil_dry * h_dry_z + gamma_soil_wet * h_wet_z)

    # Surcharge pressure
    if surcharge_type == 'Strip Load':
        depth_below_surface = max(1e-5, h_soil - z)
        res_q = stresses_stripload_retainingwall_local(
            imposedstress=q,
            width=width_surcharge,
            offset=offset_surcharge,
            toe_depth=h_soil,
            depth=depth_below_surface
        )
        sigma_q = res_q['delta sigma x [kPa]']
    else:
        sigma_q = Ka * q

    # Water pressure
    if z <= Hwtr:
        sigma_w = gamma_w * (Hwtr - z)
    else:
        sigma_w = 0.0

    # Total lateral pressure (kN/m2)
    sigma_total[node] = sigma_soil + sigma_q + sigma_w
    sigma_soil_q[node] = sigma_soil + sigma_q
    sigma_w_dict[node] = sigma_w

# Sort nodes by elevation (z) and compute tributary heights based on centroids
nodes_sorted = sorted(wall_nodes.items(), key=lambda x: x[1])
nodes_order = [n for n, _ in nodes_sorted]
z_list = [z for _, z in nodes_sorted]
Fx = {}
Fx_soil_q = {}
Fx_water_lat = {}
try:
    # build vertical boundaries: bottom = 0, internal midpoints, top = Hw
    n_stem = len(z_list)
    bounds = [0.0]
    for k in range(n_stem - 1):
        bounds.append(0.5 * (z_list[k] + z_list[k + 1]))
    bounds.append(float(Hw))

    for i, node in enumerate(nodes_order):
        h_trib = bounds[i + 1] - bounds[i]
        p_node = sigma_total.get(node, 0.0)
        p_soil_q = sigma_soil_q.get(node, 0.0)
        p_w = sigma_w_dict.get(node, 0.0)
        # lateral force on tributary area (per unit thickness t)
        Fx[node] = -p_node * h_trib * t
        Fx_soil_q[node] = -p_soil_q * h_trib * t
        Fx_water_lat[node] = -p_w * h_trib * t
except Exception:
    # fallback to previous segment approach if anything fails
    Fx = {node: 0.0 for node, _ in nodes_sorted}
    Fx_soil_q = {node: 0.0 for node, _ in nodes_sorted}
    Fx_water_lat = {node: 0.0 for node, _ in nodes_sorted}
    for i in range(len(nodes_sorted)-1):
        n1, z1 = nodes_sorted[i]
        n2, z2 = nodes_sorted[i+1]
        h = z2 - z1
        p1 = sigma_total[n1]
        p2 = sigma_total[n2]
        F = 0.5 * (p1 + p2) * h * t
        Fx[n1] += -0.5 * F
        Fx[n2] += -0.5 * F

        p1_soil_q = sigma_soil_q[n1]
        p2_soil_q = sigma_soil_q[n2]
        F_soil_q = 0.5 * (p1_soil_q + p2_soil_q) * h * t
        Fx_soil_q[n1] += -0.5 * F_soil_q
        Fx_soil_q[n2] += -0.5 * F_soil_q

        p1_w = sigma_w_dict[n1]
        p2_w = sigma_w_dict[n2]
        F_w = 0.5 * (p1_w + p2_w) * h * t
        Fx_water_lat[n1] += -0.5 * F_w
        Fx_water_lat[n2] += -0.5 * F_w

# Apply calculated earth pressure loads
for node, force in Fx.items():
    ops.load(node, force, 0.0)

#VERTICAL LOAD -------------------------------

# Apply concrete self-weight as body force to all elements
for ele in range(1, total_elements + 1):
    ops.eleLoad('-ele', ele, '-type', '-bodyForce', 0.0, -gamma_c)

#soil weight
h_dry = max(0.0, h_soil - Hwtr)
h_wet = Hwtr + h_ftg
q_soil = gamma_soil_dry * h_dry + gamma_soil_wet * h_wet   # kN/m2

# Compute tributary lengths for top-of-footing nodes (heel slab nodes)
# Use centroid-based tributary widths along the footing: domain from x=0 to x=ftg
heel_trib = {}
try:
    xs = xs_unique
    nxs = len(xs)
    # compute shared footing x-boundaries once for both top and bottom nodes
    x_bounds_footing = [0.0]
    for k in range(nxs - 1):
        x_bounds_footing.append(0.5 * (xs[k] + xs[k + 1]))
    x_bounds_footing.append(float(ftg))
    for idx, nid in enumerate(top_node_ids):
        L = x_bounds_footing[idx + 1] - x_bounds_footing[idx]
        heel_trib[nid] = L
except Exception:
    # fallback to neighbor half-distance method
    heel_trib = {}
    for idx, nid in enumerate(top_node_ids):
        if idx == 0:
            L = xs_unique[1] - xs_unique[0]
        elif idx == len(xs_unique) - 1:
            L = xs_unique[-1] - xs_unique[-2]
        else:
            L = 0.5 * (xs_unique[idx] - xs_unique[idx - 1]) + 0.5 * (xs_unique[idx + 1] - xs_unique[idx])
        heel_trib[nid] = L
# Apply soil weight and surcharge to top-of-footing nodes (heel nodes)
for idx, node in enumerate(top_node_ids):
    L = heel_trib[node]
    x_val = xs_unique[idx]
    if x_val >= (toe + bot_wall):
        Fy_soil = -q_soil * L * t
        if surcharge_type == 'Strip Load':
            d_heel = x_val - (toe + bot_wall)
            q_v = stresses_stripload_local(
                z=h_soil,
                x=d_heel - offset_surcharge,
                width=width_surcharge,
                imposedstress=q
            )['delta sigma z [kPa]']
        else:
            q_v = q
        Fy_q = -q_v * L * t
    else:
        Fy_soil = 0.0
        Fy_q = 0.0
    ops.load(node, 0.0, Fy_soil + Fy_q)
bottom_trib = {}
try:
    xs = xs_unique
    nxs = len(xs)
    x_bounds = [0.0]
    for k in range(nxs - 1):
        x_bounds.append(0.5 * (xs[k] + xs[k + 1]))
    x_bounds.append(float(ftg))
    for idx, nid in enumerate(bottom_node_ids):
        L = x_bounds_footing[idx + 1] - x_bounds_footing[idx]
        bottom_trib[nid] = L
except Exception:
    bottom_trib = {}
    for idx, nid in enumerate(bottom_node_ids):
        if idx == 0:
            L = xs_unique[1] - xs_unique[0]
        elif idx == len(xs_unique) - 1:
            L = xs_unique[-1] - xs_unique[-2]
        else:
            L = 0.5 * (xs_unique[idx] - xs_unique[idx - 1]) + 0.5 * (xs_unique[idx + 1] - xs_unique[idx])
        bottom_trib[nid] = L

# compute centroid x positions for footing nodes (centroid of each tributary interval)
centroid_x_base = {}
try:
    for idx, nid in enumerate(top_node_ids):
        cx = 0.5 * (x_bounds_footing[idx] + x_bounds_footing[idx + 1])
        centroid_x_base[nid] = cx
    for idx, nid in enumerate(bottom_node_ids):
        cx = 0.5 * (x_bounds_footing[idx] + x_bounds_footing[idx + 1])
        centroid_x_base[nid] = cx
except Exception:
    centroid_x_base = {}

uplift_pressure = {nid: gamma_w * Hwtr for nid in bottom_node_ids}

# Apply uplift to bottom nodes
for node, L in bottom_trib.items():
    Fy = uplift_pressure.get(node, 0.0) * L * t   # upward force
    ops.load(node, 0.0, Fy)

# Support information
# Fix bottom nodes in vertical DOF; fix leftmost bottom node in both DOFs
leftmost_bottom = bottom_node_ids[0]
ops.fix(leftmost_bottom, 1, 1)
for nid in bottom_node_ids[1:]:
    ops.fix(nid, 0, 1)
# Also fix the top-of-footing leftmost node in horizontal DOF to prevent rigid
ops.fix(top_node_ids[0], 1, 0)

ops.system('BandSPD')
ops.numberer('RCM')
ops.constraints('Plain')
ops.integrator('LoadControl', 1.0)

# Analysis type
ops.analysis('Static')
ops.analyze(1) # Specifies how many times the analysis will be performed.

# ---------------------------------
# Output
# ---------------------------------

# --- Seismic treatment ---
# Option A: Pseudo-static (apply kh*W lateral load)
# --- Seismic treatment (Pseudo-static automatically included) ---
lateral_per_node = 0.0
kh = PGA * FPGA
if kh > 0.0:
    st.info('Applying pseudo-static seismic loads...')
    # approximate weights (kN) per unit thickness * thickness `t`
    wall_area = (top_wall + bot_wall) / 2.0 * Hw
    footing_area = ftg * h_ftg
    W_conc = gamma_c * (wall_area + footing_area) * t
    h_dry = max(0.0, h_soil - Hwtr)
    h_wet = Hwtr + h_ftg
    W_soil = (gamma_soil_dry * h_dry + gamma_soil_wet * h_wet) * ftg * t
    W_total = W_soil + W_conc

    total_lateral = kh * W_total
    # distribute equally to back-face wall nodes
    nodes = list(wall_nodes.keys())
    n_w = len(nodes) if len(nodes) > 0 else 1
    lateral_per_node = - total_lateral / n_w

    # Apply as a new plain load pattern so existing loads remain
    pseudo_pattern = 3
    ops.pattern('Plain', pseudo_pattern, 1)
    for node in nodes:
        ops.load(node, lateral_per_node, 0.0)

    # Re-run a static analysis to obtain displacements under combined loads
    ops.system('BandSPD')
    ops.numberer('RCM')
    ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.analysis('Static')
    ok_ps = ops.analyze(1)
    if ok_ps == 0:
        max_x = 0.0
        max_y = 0.0
        for i in range(1, total_nodes + 1):
            try:
                dx = abs(ops.nodeDisp(i, 1))
                dy = abs(ops.nodeDisp(i, 2))
                if dx > max_x:
                    max_x = dx
                if dy > max_y:
                    max_y = dy
            except Exception:
                pass
        st.success(f'Pseudo-static applied (kh={kh:.4f}). Total lateral = {total_lateral:.3f} kN. Max x disp = {max_x:.6f} m')
    else:
        st.error('Pseudo-static analysis failed (ops.analyze returned non-zero).')

# Displacements
for i in range(1, total_nodes + 1):
    try:
        print(i, "node x=", ops.nodeDisp(i, 1), "m  -  y=", ops.nodeDisp(i, 2), "m")
    except Exception:
        pass


# =========================================================
# CONSOLIDATED COMPUTATIONS (Stability, Loads, and Metrics)
# =========================================================

# 1. Total lateral (resultant) from Fx dictionary (kN)
try:
    total_lateral_resultant = -sum(Fx.values())  # sign: Fx stored as negative for earth pressure
except Exception:
    total_lateral_resultant = 0.0

try:
    soil_node_loads = {}
    for idx, nid in enumerate(top_node_ids):
        L = heel_trib[nid]
        x_val = xs_unique[idx]
        if x_val >= (toe + bot_wall):
            if surcharge_type == 'Strip Load':
                d_heel = x_val - (toe + bot_wall)
                q_v = stresses_stripload_local(
                    z=h_soil,
                    x=d_heel - offset_surcharge,
                    width=width_surcharge,
                    imposedstress=q
                )['delta sigma z [kPa]']
            else:
                q_v = q
            soil_node_loads[nid] = -(q_soil + q_v) * L * t
        else:
            soil_node_loads[nid] = 0.0
    soil_lines = "\n".join([f"- Node {n}: {soil_node_loads[n]:.3f} kN" for n in sorted(soil_node_loads.keys())])
    soil_total = sum(soil_node_loads.values())
except Exception:
    soil_lines = "(not available)"
    soil_total = 0.0
    soil_node_loads = {}

# 3. Uplift applied at bottom nodes (bottom_trib)
try:
    uplift_node_loads = {n: (uplift_pressure.get(n, 0.0) * L * t) for n, L in bottom_trib.items()}
    uplift_lines = "\n".join([f"- Node {n}: {uplift_node_loads[n]:.3f} kN" for n in sorted(uplift_node_loads.keys())])
    uplift_total = sum(uplift_node_loads.values())
except Exception:
    uplift_lines = "(not available)"
    uplift_total = 0.0
    uplift_node_loads = {}

net_vertical = soil_total + uplift_total

# 4. Concrete weights and lever arms
base_area = ftg * h_ftg
W_base = gamma_c * base_area * t
cg_x_base = ftg / 2.0

stem_rect_area = top_wall * Hw
W_stem_rect = gamma_c * stem_rect_area * t
cg_x_stem_rect = toe + bot_wall - top_wall / 2.0

stem_tri_area = 0.5 * taper * Hw
W_stem_tri = gamma_c * stem_tri_area * t
cg_x_stem_tri = toe + (2.0 / 3.0) * taper if taper > 0 else float(toe)

W_conc = W_base + W_stem_rect + W_stem_tri
wall_area = stem_rect_area + stem_tri_area
footing_area = base_area

# 5. Soil weight on heel and lever arm
soil_width = ftg - (toe + bot_wall)
W_soil = (gamma_soil_dry * h_dry + gamma_soil_wet * h_wet) * soil_width * t
cg_x_soil = toe + bot_wall + soil_width / 2.0

# 5.5 Surcharge vertical weight and lever arm on heel footing
if surcharge_type == 'Strip Load':
    n_pts = 10
    forces_pts = []
    for k in range(n_pts):
        frac = (k + 0.5) / n_pts
        d_pt = frac * soil_width
        x_rel = d_pt - offset_surcharge
        res = stresses_stripload_local(z=h_soil, x=x_rel, width=width_surcharge, imposedstress=q)
        forces_pts.append(max(0.0, res['delta sigma z [kPa]']))
    avg_q_v = np.mean(forces_pts)
else:
    avg_q_v = q
W_surcharge = avg_q_v * soil_width * t

# 6. Global stability checks (about x_pivot = 0.0)
x_pivot_global = 0.0
try:
    M_ot_global = 0.0
    for n, f in Fx.items():
        try:
            coord = ops.nodeCoord(n)
            z_n = float(coord[1])
        except Exception:
            z_n = 0.0
        M_ot_global += abs(f) * abs(z_n)
except Exception:
    M_ot_global = 0.0

try:
    M_uplift_global = uplift_total * (ftg / 2.0 - x_pivot_global)
except Exception:
    M_uplift_global = 0.0

M_res_global = (W_base * abs(cg_x_base - x_pivot_global) +
                W_stem_rect * abs(cg_x_stem_rect - x_pivot_global) +
                W_stem_tri * abs(cg_x_stem_tri - x_pivot_global) +
                (W_soil + W_surcharge) * abs(cg_x_soil - x_pivot_global) -
                M_uplift_global)

try:
    W_down_global = W_conc + W_soil + W_surcharge - uplift_total
except Exception:
    W_down_global = W_conc + W_soil + W_surcharge

try:
    FS_ot_global = M_res_global / M_ot_global if M_ot_global > 0 else float('inf')
except Exception:
    FS_ot_global = None

try:
    phi_rad = np.radians(phi)
    R_slide = max(0.0, W_down_global) * np.tan(phi_rad)
    F_drive = total_lateral_resultant
    FS_slide_global = R_slide / F_drive if F_drive > 0 else float('inf')
except Exception:
    R_slide = 0.0
    FS_slide_global = None

# 7. Stem stability checks (about x_pivot = toe)
x_pivot_stem = float(toe)
x_centroid = float(ftg) / 2.0
try:
    M_ot_stem = 0.0
    for n, f in Fx.items():
        if n in centroid_x_wall:
            x_n = centroid_x_wall[n]
        elif n in centroid_x_base:
            x_n = centroid_x_base[n]
        else:
            try:
                coord = ops.nodeCoord(n)
                x_n = float(coord[0])
            except Exception:
                x_n = x_centroid
        M_ot_stem += (abs(f)) * abs(x_n - x_pivot_stem)
except Exception:
    M_ot_stem = 0.0

try:
    W_down_stem = W_conc + W_soil - uplift_total
except Exception:
    W_down_stem = W_conc + W_soil
M_res_stem = W_down_stem * abs(x_centroid - x_pivot_stem)
try:
    FS_ot_stem = M_res_stem / M_ot_stem if M_ot_stem > 0 else float('inf')
except Exception:
    FS_ot_stem = None

# 8. Seismic coefficients and peak estimates
kh = PGA * FPGA
F_seismic_pseudo = kh * (W_soil + W_conc)
g_acc = 9.81
peak_accel = PGA * g_acc * FPGA
inertial_peak = (W_soil + W_conc) * (peak_accel / g_acc)

# 9. Get maximum displacements
max_disp_x = 0.0
max_disp_y = 0.0
for i in range(1, total_nodes + 1):
    try:
        dx = abs(ops.nodeDisp(i, 1))
        dy = abs(ops.nodeDisp(i, 2))
        if dx > max_disp_x:
            max_disp_x = dx
        if dy > max_disp_y:
            max_disp_y = dy
    except Exception:
        pass

# 10. Generate Summary DataFrames
df_lateral = pd.DataFrame()
df_vertical = pd.DataFrame()
try:
    # Lateral table (wall nodes)
    lateral_rows = []
    total_lat_force = 0.0
    total_lat_moment = 0.0
    for n in sorted(Fx.keys()):
        f = Fx[n]
        try:
            coord = ops.nodeCoord(n)
            z_n = float(coord[1])
        except Exception:
            z_n = 0.0
        moment = abs(f) * abs(z_n)
        lateral_rows.append({'Node': n, 'Y (m) lever': z_n, 'Force (kN)': f, 'Moment (kN·m)': moment})
        total_lat_force += f
        total_lat_moment += moment

    df_lateral = pd.DataFrame(lateral_rows)
    df_lateral = pd.concat([df_lateral, pd.DataFrame([{'Node': 'TOTAL', 'Y (m) lever': np.nan, 'Force (kN)': total_lat_force, 'Moment (kN·m)': total_lat_moment}])], ignore_index=True)
    df_lateral['Node'] = df_lateral['Node'].astype(str)

    # Vertical table (top nodes soil and bottom uplift)
    vertical_rows = []
    total_vert_force = 0.0
    total_vert_moment = 0.0
    vert_nodes = set(list(heel_trib.keys()) + list(bottom_trib.keys()))
    for n in sorted(vert_nodes):
        soil_f = abs(soil_node_loads.get(n, 0.0))
        uplift_f = -abs(uplift_node_loads.get(n, 0.0))
        f_vert = soil_f + uplift_f
        if n in centroid_x_base:
            x_n = centroid_x_base[n]
        else:
            try:
                coord = ops.nodeCoord(n)
                x_n = float(coord[0])
            except Exception:
                x_n = float(ftg) / 2.0
        moment = f_vert * abs(x_n - x_pivot_global)
        vertical_rows.append({'Node': n, 'X (m)': x_n, 'Soil (kN)': soil_node_loads.get(n, 0.0), 'Uplift (kN)': uplift_node_loads.get(n, 0.0), 'Total Vertical (kN)': f_vert, 'Moment (kN·m)': moment})
        total_vert_force += f_vert
        total_vert_moment += moment

    df_vertical = pd.DataFrame(vertical_rows)
    df_vertical = pd.concat([df_vertical, pd.DataFrame([{'Node': 'TOTAL', 'X (m)': np.nan, 'Soil (kN)': np.nan, 'Uplift (kN)': np.nan, 'Total Vertical (kN)': total_vert_force, 'Moment (kN·m)': total_vert_moment}])], ignore_index=True)
    df_vertical['Node'] = df_vertical['Node'].astype(str)
except Exception as e:
    pass

# Helper to get stem moments
def get_stem_moments():
    node_moments = {}
    
    # 1. Base of the stem
    base_nodes = [left_node, right_node]
    M_base = 0.0
    centroid_base = (ops.nodeCoord(left_node, 1) + ops.nodeCoord(right_node, 1))/2.0
    for ele in stem_elements_by_layer[0]:
        nodes = ops.eleNodes(ele)
        forces = ops.eleResponse(ele, 'forces')
        for i, node in enumerate(nodes):
            if node in base_nodes:
                Fy = forces[i*2 + 1]
                x = ops.nodeCoord(node, 1)
                M_base -= Fy * (x - centroid_base)
    for node in base_nodes:
        node_moments[node] = M_base
    
    # 2. Layer cuts (1 to 8)
    for j in range(1, 9):
        b_j, f_j = stem_pairs[j-1]
        cut_nodes = [b_j, f_j]
        M_cut = 0.0
        centroid_cut = (ops.nodeCoord(b_j, 1) + ops.nodeCoord(f_j, 1))/2.0
        for ele in stem_elements_by_layer[j-1]:
            nodes = ops.eleNodes(ele)
            forces = ops.eleResponse(ele, 'forces')
            for i, node in enumerate(nodes):
                if node in cut_nodes:
                    Fy = forces[i*2 + 1]
                    x = ops.nodeCoord(node, 1)
                    M_cut += Fy * (x - centroid_cut)
        for node in cut_nodes:
            node_moments[node] = M_cut
        
    return node_moments


# =========================================================
# CONCRETE PROPERTIES ULTIMATE BENDING CAPACITY CODE
# =========================================================

# Concrete Stress Block parameters (AS3600 approximation)
alpha_stress = max(0.67, min(0.85, 0.85 - 0.0015 * fc))
gamma_stress = max(0.67, min(0.85, 0.97 - 0.0025 * fc))

try:
    from concreteproperties.material import Concrete, SteelBar
    import concreteproperties.stress_strain_profile as ssp
    from concreteproperties.pre import add_bar_rectangular_array
    from concreteproperties.concrete_section import ConcreteSection
    from sectionproperties.pre.library.primitive_sections import rectangular_section
    
    # Define Materials
    material_concrete = Concrete(
        name="Concrete Class",
        density=2.4e-6,
        stress_strain_profile=ssp.ConcreteLinear(elastic_modulus=Ec),
        ultimate_stress_strain_profile=ssp.RectangularStressBlock(
            compressive_strength=fc,
            alpha=alpha_stress,
            gamma=gamma_stress,
            ultimate_strain=0.003,
        ),
        flexural_tensile_strength=ft_tensile,
        colour="lightgrey",
    )
    
    material_steel = SteelBar(
        name="Steel Class",
        density=7.85e-6,
        stress_strain_profile=ssp.SteelElasticPlastic(
            yield_strength=fy,
            elastic_modulus=Es,
            fracture_strain=0.05,
        ),
        colour="grey",
    )
    
    # 1. Stem Section Bending Capacity (1m wide strip)
    stem_beam_width = 1000.0  # mm
    stem_beam_height = bot_wall * 1000.0  # mm
    
    geom_stem = rectangular_section(d=stem_beam_height, b=stem_beam_width, material=material_concrete)
    n_x_bars_stem = int((stem_beam_width - 2 * stem_cover) / stem_spacing_x) + 1
    
    # Bottom bars
    geom_stem = add_bar_rectangular_array(
        geometry=geom_stem,
        area=stem_bar_area,
        material=material_steel,
        n_x=n_x_bars_stem,
        x_s=stem_spacing_x,
        anchor=(stem_cover, stem_cover),
    )
    
    # Top bars
    geom_stem = add_bar_rectangular_array(
        geometry=geom_stem,
        area=stem_bar_area,
        material=material_steel,
        n_x=n_x_bars_stem,
        x_s=stem_spacing_x,
        anchor=(stem_cover, stem_beam_height - stem_cover),
    )
    
    # Side bars
    geom_stem = add_bar_rectangular_array(
        geometry=geom_stem,
        area=stem_sidebar_area,
        material=material_steel,
        n_x=2,
        x_s=stem_beam_width - 2 * stem_cover,
        n_y=3,
        y_s=stem_spacing_y,
        anchor=(stem_cover, stem_cover + 150),
    )
    
    conc_sec_stem = ConcreteSection(geom_stem)
    stem_capacity_sag = conc_sec_stem.ultimate_bending_capacity(theta=0)
    stem_capacity_hog = conc_sec_stem.ultimate_bending_capacity(theta=np.pi)
    
    M_u_stem_allow = min(stem_capacity_sag.m_xy, stem_capacity_hog.m_xy) / 1e6  # kN.m
    
    # Generate Stem section geometry plot
    plt.close('all')
    conc_sec_stem.plot_section()
    fig_stem_sec = plt.gcf()
    fig_stem_sec.set_size_inches(6, 4)
    # Set labels
    try:
        ax = fig_stem_sec.axes[0]
        ax.set_title("Stem Reinforcement Cross Section (1m Width)", fontsize=10)
    except Exception:
        pass
    fig_stem_sec.tight_layout()
    plt.close(fig_stem_sec)

    # 2. Footing Section Bending Capacity (1m wide strip)
    ftg_beam_width = 1000.0  # mm
    ftg_beam_height = h_ftg * 1000.0  # mm
    
    geom_ftg = rectangular_section(d=ftg_beam_height, b=ftg_beam_width, material=material_concrete)
    n_x_bars_ftg = int((ftg_beam_width - 2 * ftg_cover) / ftg_spacing_x) + 1
    n_y_bars_ftg = int((ftg_beam_height - 2 * ftg_cover) / ftg_spacing_y) + 1
    
    # Bottom bars
    geom_ftg = add_bar_rectangular_array(
        geometry=geom_ftg,
        area=ftg_bar_area,
        material=material_steel,
        n_x=n_x_bars_ftg,
        x_s=ftg_spacing_x,
        anchor=(ftg_cover, ftg_cover),
    )
    
    # Top bars
    geom_ftg = add_bar_rectangular_array(
        geometry=geom_ftg,
        area=ftg_bar_area,
        material=material_steel,
        n_x=n_x_bars_ftg,
        x_s=ftg_spacing_x,
        anchor=(ftg_cover, ftg_beam_height - ftg_cover),
    )
    
    # Side bars
    geom_ftg = add_bar_rectangular_array(
        geometry=geom_ftg,
        area=ftg_sidebar_area,
        material=material_steel,
        n_x=2,
        x_s=ftg_beam_width - 2 * ftg_cover,
        n_y=n_y_bars_ftg,
        y_s=ftg_spacing_y,
        anchor=(ftg_cover, ftg_cover),
    )
    
    conc_sec_ftg = ConcreteSection(geom_ftg)
    ftg_capacity_sag = conc_sec_ftg.ultimate_bending_capacity(theta=0)
    ftg_capacity_hog = conc_sec_ftg.ultimate_bending_capacity(theta=np.pi)
    
    M_u_ftg_sag_allow = ftg_capacity_sag.m_xy / 1e6  # kN.m (sagging, tension at bottom)
    M_u_ftg_hog_allow = ftg_capacity_hog.m_xy / 1e6  # kN.m (hogging, tension at top)
    
    # Generate Footing section geometry plot
    plt.close('all')
    conc_sec_ftg.plot_section()
    fig_ftg_sec = plt.gcf()
    fig_ftg_sec.set_size_inches(6, 4)
    try:
        ax = fig_ftg_sec.axes[0]
        ax.set_title("Footing Reinforcement Cross Section (1m Width)", fontsize=10)
    except Exception:
        pass
    fig_ftg_sec.tight_layout()
    plt.close(fig_ftg_sec)

except Exception as e:
    st.sidebar.error(f"Error calculating bending capacity: {e}")
    M_u_stem_allow = 0.0
    M_u_ftg_sag_allow = 0.0
    M_u_ftg_hog_allow = 0.0
    fig_stem_sec = None
    fig_ftg_sec = None

# ----------------------------------------------------
# Actual Moment Demands Calculations
# ----------------------------------------------------
try:
    # Stem Wall actual moment from FEM base moment
    node_moments_dict = get_stem_moments()
    M_actual_stem = abs(node_moments_dict[left_node]) / t
except Exception:
    M_actual_stem = 0.0

try:
    B_len = ftg
    Rv_force = max(0.1, W_conc + W_soil + W_surcharge - uplift_total)
    M_net_moment = M_res_global - M_ot_global
    x_R_loc = M_net_moment / Rv_force
    e_ecc = B_len / 2.0 - x_R_loc

    if abs(e_ecc) <= B_len / 6.0:
        q_t = (Rv_force / B_len) * (1.0 + 6.0 * e_ecc / B_len)
        q_h = (Rv_force / B_len) * (1.0 - 6.0 * e_ecc / B_len)
    elif e_ecc > B_len / 6.0:
        q_t = (2.0 * Rv_force) / (3.0 * x_R_loc)
        q_h = 0.0
    else:
        q_h = (2.0 * Rv_force) / (3.0 * (B_len - x_R_loc))
        q_t = 0.0

    def get_q_bearing_at(x):
        if abs(e_ecc) <= B_len / 6.0:
            return q_t + (q_h - q_t) * (x / B_len)
        elif e_ecc > B_len / 6.0:
            L_b = 3.0 * x_R_loc
            return q_t * (1.0 - x / L_b) if x <= L_b else 0.0
        else:
            L_b = 3.0 * (B_len - x_R_loc)
            return q_h * (x - (B_len - L_b)) / L_b if x >= B_len - L_b else 0.0

    q_toe_end = get_q_bearing_at(toe)
    q_heel_start = get_q_bearing_at(toe + bot_wall)

    # Toe actual moment (tension at bottom, sagging)
    w_down_t = (gamma_c * h_ftg + gamma_soil_dry * h_soil_toe) * t
    M_bearing_t = (get_q_bearing_at(0.0) * (toe**2) / 3.0 + q_toe_end * (toe**2) / 6.0) * t
    M_down_t = w_down_t * (toe**2) / 2.0
    M_actual_toe = max(0.0, M_bearing_t - M_down_t) / t

    # Heel actual moment (tension at top, hogging)
    w_down_h = (gamma_c * h_ftg + (gamma_soil_dry * h_dry + gamma_soil_wet * h_wet) / B_len * soil_width + q) * t
    heel_len = B_len - (toe + bot_wall)
    M_bearing_h = (q_heel_start * (heel_len**2) / 6.0 + get_q_bearing_at(B_len) * (heel_len**2) / 3.0) * t
    M_down_h = w_down_h * (heel_len**2) / 2.0
    M_actual_heel = max(0.0, M_down_h - M_bearing_h) / t
    sigma_max = max(q_t, q_h)
except Exception:
    M_actual_toe = 0.0
    M_actual_heel = 0.0
    sigma_max = 0.0

# Compute overburden stresses and shear force at footing base (used by bearing and pile analysis)
D_depth = h_ftg + h_soil_toe
H_water_front = h_ftg + Hwtr_front if Hwtr_front > 0 else 0.0
if H_water_front >= D_depth:
    p0_eff = D_depth * (gamma_soil_wet - gamma_w)
    base_sigma_v = D_depth * gamma_soil_wet
elif H_water_front > 0:
    p0_eff = H_water_front * (gamma_soil_wet - gamma_w) + (D_depth - H_water_front) * gamma_soil_dry
    base_sigma_v = H_water_front * gamma_soil_wet + (D_depth - H_water_front) * gamma_soil_dry
else:
    p0_eff = D_depth * gamma_soil_dry
    base_sigma_v = D_depth * gamma_soil_dry

p0_eff = max(0.0, p0_eff)
base_sigma_v = max(0.0, base_sigma_v)
H_force = abs(F_drive) if 'F_drive' in locals() else 0.0

# ----------------------------------------------------
# Bearing Capacity Analysis using Groundhog
# ----------------------------------------------------
q_ult = 0.0
q_allow = 0.0
FS_actual_bearing = float('inf')
bearing_pass = True

if enable_bearing:
    try:
        import groundhog.shallowfoundations.capacity as cap
        
        # Calculate B' and parameters for groundhog
        B_prime = max(0.01, ftg - 2.0 * abs(e_ecc)) if 'e_ecc' in locals() else ftg
        L_prime = 1e5

        # Unit weight below foundation base
        if Hwtr_front > 0 or Hwtr > 0:
            effective_unit_weight_below = max(3.0, gamma_soil_wet - gamma_w)
        else:
            effective_unit_weight_below = gamma_soil_dry

        # Horizontal and vertical forces
        V_force = max(0.1, W_down_global) if 'W_down_global' in locals() else 0.1

        if bearing_soil_type == 'Sand (Drained)':
            load_inc = np.degrees(np.arctan(H_force / V_force))
            # Call groundhog drained capacity
            res = cap.verticalcapacity_drained_api(
                vertical_effective_stress=p0_eff,
                effective_friction_angle=float(phi),
                effective_unit_weight=effective_unit_weight_below,
                effective_length=L_prime,
                effective_width=B_prime,
                base_depth=D_depth,
                skirted=False,
                load_inclination=float(load_inc),
                validate=False
            )
            q_ult = res.get('qu [kPa]', 0.0)
            if np.isnan(q_ult):
                q_ult = 0.0
        else:
            # Clay (Undrained)
            res = cap.verticalcapacity_undrained_api(
                effective_length=L_prime,
                effective_width=B_prime,
                su_base=su_val,
                base_depth=D_depth,
                skirted=False,
                base_sigma_v=base_sigma_v,
                horizontal_load=H_force * L_prime,
                validate=False
            )
            q_ult = res.get('qu [kPa]', 0.0)
            if np.isnan(q_ult):
                q_ult = 0.0

        q_allow = q_ult / FS_bearing
        if sigma_max > 0:
            FS_actual_bearing = q_ult / sigma_max
        else:
            FS_actual_bearing = float('inf')
            
        bearing_pass = FS_actual_bearing >= FS_bearing
    except Exception as e:
        q_ult = 0.0
        q_allow = 0.0
        FS_actual_bearing = 0.0
        bearing_pass = False



# ----------------------------------------------------
# Base Anchor (Pile) Foundation Analysis
# ----------------------------------------------------
pile_pass_toe_axial = True
pile_pass_toe_lateral = True
pile_pass_heel_axial = True
pile_pass_heel_lateral = True
pile_pass_toe_interaction = True
pile_pass_heel_interaction = True
fig_pile_sec = None

if enable_pile:
    try:
        from openpile.construct import PileSection, Pile, SoilProfile, Layer, Model
        from openpile.soilmodels import API_sand, API_clay
        from groundhog.deepfoundations.axialcapacity.skinfriction import API_unit_shaft_friction_sand_rp2geo, API_unit_shaft_friction_clay
        from groundhog.deepfoundations.axialcapacity.endbearing import API_unit_end_bearing_sand_rp2geo, API_unit_end_bearing_clay
        from sectionproperties.pre.geometry import CompoundGeometry
        from sectionproperties.pre.library import circular_section, rectangular_section
        from concreteproperties import (
            Concrete,
            ConcreteLinear,
            ConcreteSection,
            RectangularStressBlock,
            SteelBar,
            SteelElasticPlastic,
        )
        from concreteproperties.pre import add_bar_circular_array, add_bar_rectangular_array
        import matplotlib.path as mpath
        
        # 1. Custom Rectangular Section for OpenPile
        class RectangularPileSection(PileSection):
            top: float
            bottom: float
            width_x: float
            width_y: float

            @property
            def top_elevation(self) -> float:
                return self.top

            @property
            def bottom_elevation(self) -> float:
                return self.bottom

            @property
            def length(self) -> float:
                return self.top - self.bottom

            @property
            def area(self) -> float:
                return self.width_x * self.width_y

            @property
            def footprint(self) -> float:
                return self.width_x * self.width_y

            @property
            def outer_perimeter(self) -> float:
                return 2.0 * (self.width_x + self.width_y)

            @property
            def inner_perimeter(self) -> float:
                return 0.0

            @property
            def entrapped_area(self) -> float:
                return 0.0

            @property
            def second_moment_of_area(self) -> float:
                return self.width_y * (self.width_x ** 3) / 12.0

            @property
            def width(self) -> float:
                return self.width_x

            def get_volume(self, length) -> float:
                return self.area * length

            def get_entrapped_volume(self, length) -> float:
                return 0.0

        # 2. Force split calculations
        d_piles = ftg - 2.0 * pile_offset
        P_toe = Rv_force * (ftg - pile_offset - x_R_loc) / d_piles
        P_heel = Rv_force * (x_R_loc - pile_offset) / d_piles
        H_toe = H_force / 2.0
        H_heel = H_force / 2.0

        # 3. Single pile capacities (using groundhog)
        gamma_sub = gamma_pile_soil - gamma_w if (Hwtr_front > 0 or Hwtr > 0) else gamma_pile_soil
        p0_mid = p0_eff + (L_pile / 2.0) * gamma_sub
        p0_tip = p0_eff + L_pile * gamma_sub

        if pile_soil_type == 'Sand (Drained)':
            if phi_pile_soil < 36:
                rd_cat = "Medium dense"
            elif phi_pile_soil < 40:
                rd_cat = "Dense"
            else:
                rd_cat = "Very dense"
            sf_res = API_unit_shaft_friction_sand_rp2geo(api_relativedensity=rd_cat, api_soildescription='Sand', sigma_vo_eff=p0_mid, validate=False)
            f_s_comp = sf_res['f_s_comp_out [kPa]']
            f_s_tens = sf_res['f_s_tens_out [kPa]']
            eb_res = API_unit_end_bearing_sand_rp2geo(api_relativedensity=rd_cat, api_soildescription='Sand', sigma_vo_eff=p0_tip, validate=False)
            q_b = eb_res['q_b_plugged [kPa]']
        else:
            sf_res = API_unit_shaft_friction_clay(undrained_shear_strength=su_pile_soil, sigma_vo_eff=p0_mid, validate=False)
            f_s_comp = sf_res['f_s_comp_out [kPa]']
            f_s_tens = sf_res['f_s_tens_out [kPa]']
            eb_res = API_unit_end_bearing_clay(undrained_shear_strength=su_pile_soil, validate=False)
            q_b = eb_res['q_b_plugged [kPa]']

        if pile_shape == 'Circle':
            A_base = np.pi * (diameter_pile ** 2) / 4.0
            A_shaft = np.pi * diameter_pile * L_pile
        else:
            A_base = width_x_pile * width_y_pile
            A_shaft = 2.0 * (width_x_pile + width_y_pile) * L_pile

        Q_shaft_comp = f_s_comp * A_shaft
        Q_shaft_tens = f_s_tens * A_shaft
        Q_base = q_b * A_base
        Q_comp_ult = Q_shaft_comp + Q_base
        Q_tens_ult = Q_shaft_tens
        Q_comp_allow = Q_comp_ult / FS_pile_axial
        Q_tens_allow = Q_tens_ult / FS_pile_axial

        # Broms' Lateral Capacity
        if pile_soil_type == 'Sand (Drained)':
            Kp = np.tan(np.radians(45.0 + phi_pile_soil / 2.0)) ** 2
            H_ult = 0.5 * gamma_pile_soil * (diameter_pile if pile_shape == 'Circle' else width_x_pile) * (L_pile ** 3) * Kp
        else:
            H_ult = 9.0 * su_pile_soil * (diameter_pile if pile_shape == 'Circle' else width_x_pile) * (L_pile - 1.5 * (diameter_pile if pile_shape == 'Circle' else width_x_pile))
        H_allow = H_ult / FS_pile_lateral

        # 4. OpenPile Lateral Winkler Analysis
        if pile_shape == 'Circle':
            p_elem = Pile.create_tubular(
                name="Circular Pile",
                top_elevation=0.0,
                bottom_elevation=-L_pile,
                diameter=diameter_pile,
                wt=diameter_pile / 2.0
            )
        else:
            p_elem = Pile(
                name="Rectangular Pile",
                material="Concrete",
                sections=[
                    RectangularPileSection(
                        top=0.0,
                        bottom=-L_pile,
                        width_x=width_x_pile,
                        width_y=width_y_pile
                    )
                ]
            )

        if pile_soil_type == 'Sand (Drained)':
            lat_model = API_sand(phi=phi_pile_soil, kind=pile_loading_type)
        else:
            lat_model = API_clay(Su=[su_pile_soil, su_pile_soil], eps50=eps50_pile_soil, kind=pile_loading_type)

        sp_elem = SoilProfile(
            name="Pile Soil Profile",
            top_elevation=0.0,
            water_line=0.0 if (Hwtr_front > 0 or Hwtr > 0) else -L_pile,
            layers=[
                Layer(
                    name="layer1",
                    top=0.0,
                    bottom=-L_pile,
                    weight=gamma_pile_soil,
                    lateral_model=lat_model
                )
            ]
        )

        # Winkler solve for toe
        M_toe_model = Model(name="Toe Pile Model", pile=p_elem, soil=sp_elem, element_type="EulerBernoulli")
        M_toe_model.set_support(elevation=-L_pile, Tz=True)
        M_toe_model.set_pointload(elevation=0.0, Py=H_toe, Pz=-P_toe)
        res_toe = M_toe_model.solve()

        # Winkler solve for heel
        M_heel_model = Model(name="Heel Pile Model", pile=p_elem, soil=sp_elem, element_type="EulerBernoulli")
        M_heel_model.set_support(elevation=-L_pile, Tz=True)
        M_heel_model.set_pointload(elevation=0.0, Py=H_heel, Pz=-P_heel)
        res_heel = M_heel_model.solve()

        # Extract winkler results
        V_max_toe = max(abs(res_toe.forces['V [kN]']))
        M_max_toe = max(abs(res_toe.forces['M [kNm]']))
        V_max_heel = max(abs(res_heel.forces['V [kN]']))
        M_max_heel = max(abs(res_heel.forces['M [kNm]']))

        # 5. concrete-properties Interaction Diagram
        material_concrete_pile = Concrete(
            name="Pile Concrete",
            density=2.4e-6,
            stress_strain_profile=ConcreteLinear(elastic_modulus=Ec),
            ultimate_stress_strain_profile=RectangularStressBlock(
                compressive_strength=fc,
                alpha=0.85,
                gamma=0.85,
                ultimate_strain=0.003
            ),
            flexural_tensile_strength=3.0,
            colour="lightgrey"
        )
        material_steel_pile = SteelBar(
            name="Pile Rebar",
            density=7.85e-6,
            stress_strain_profile=SteelElasticPlastic(
                yield_strength=fy,
                elastic_modulus=Es,
                fracture_strain=0.05
            ),
            colour="grey"
        )

        if pile_material == 'Concrete':
            if pile_shape == 'Circle':
                geom_pile = circular_section(d=diameter_pile * 1000.0, n=32, material=material_concrete_pile)
                rebar_area = np.pi * (rebar_dia_pile ** 2) / 4.0
                geom_pile = add_bar_circular_array(
                    geometry=geom_pile,
                    area=rebar_area,
                    material=material_steel_pile,
                    n_bar=int(n_rebar_pile),
                    r_array=diameter_pile * 1000.0 / 2.0 - cover_pile,
                    ctr=(0.0, 0.0)
                )
            else:
                geom_pile = rectangular_section(d=width_x_pile * 1000.0, b=width_y_pile * 1000.0, material=material_concrete_pile)
                rebar_area = np.pi * (rebar_dia_pile ** 2) / 4.0
                geom_pile = add_bar_rectangular_array(
                    geometry=geom_pile,
                    area=rebar_area,
                    material=material_steel_pile,
                    n_x=2,
                    x_s=width_y_pile * 1000.0 - 2.0 * cover_pile,
                    n_y=2,
                    y_s=width_x_pile * 1000.0 - 2.0 * cover_pile,
                    anchor=(cover_pile, cover_pile)
                )
        else:
            # Steel solid section model
            steel_concrete = Concrete(
                name="Solid Steel",
                density=7.85e-6,
                stress_strain_profile=ConcreteLinear(elastic_modulus=Es),
                ultimate_stress_strain_profile=RectangularStressBlock(
                    compressive_strength=fy,
                    alpha=1.0,
                    gamma=1.0,
                    ultimate_strain=0.05
                ),
                flexural_tensile_strength=fy,
                colour="grey"
            )
            if pile_shape == 'Circle':
                geom_pile = circular_section(d=diameter_pile * 1000.0, n=32, material=steel_concrete)
            else:
                geom_pile = rectangular_section(d=width_x_pile * 1000.0, b=width_y_pile * 1000.0, material=steel_concrete)

        if not isinstance(geom_pile, CompoundGeometry):
            geom_pile = CompoundGeometry([geom_pile])
        sec_pile = ConcreteSection(geom_pile)

        # Generate Pile section geometry plot
        plt.close('all')
        sec_pile.plot_section()
        fig_pile_sec = plt.gcf()
        fig_pile_sec.set_size_inches(5, 4)
        try:
            ax = fig_pile_sec.axes[0]
            ax.set_title("Pile Cross Section Geometry", fontsize=10)
        except Exception:
            pass
        fig_pile_sec.tight_layout()
        plt.close(fig_pile_sec)

        mi_res = sec_pile.moment_interaction_diagram(control_points=[('kappa0', 0.0), ('d_n', 1e-6)], progress_bar=False)
        pos_m = [r.m_x / 1.0e6 for r in mi_res.results]
        pos_n = [r.n / 1000.0 for r in mi_res.results]

        curve_m = pos_m + [-x for x in pos_m[::-1]]
        curve_n = pos_n + pos_n[::-1]
        polygon_points = list(zip(curve_m, curve_n))
        path_pile = mpath.Path(polygon_points)

        # Verification Checks
        pile_pass_toe_axial = P_toe <= Q_comp_allow if P_toe >= 0 else abs(P_toe) <= Q_tens_allow
        pile_pass_toe_lateral = V_max_toe <= H_allow
        pile_pass_toe_interaction = path_pile.contains_point((M_max_toe, P_toe))

        pile_pass_heel_axial = P_heel <= Q_comp_allow if P_heel >= 0 else abs(P_heel) <= Q_tens_allow
        pile_pass_heel_lateral = V_max_heel <= H_allow
        pile_pass_heel_interaction = path_pile.contains_point((M_max_heel, P_heel))

    except Exception as e:
        st.error(f"Error during Pile Analysis: {e}")
        import traceback
        st.code(traceback.format_exc())

# =========================================================
# DASHBOARD UI RENDER
# =========================================================

# 1. Title Banner
st.markdown("""
<div class="dashboard-header">
    <h1>🧱 Cantilever Retaining Wall Designer</h1>
    <p>Preliminary Cantilever Retaining Wall Design tool</p>
</div>
""", unsafe_allow_html=True)

# 2. KPI Summary Row
if enable_bearing:
    kpi_cols = st.columns(6)
else:
    kpi_cols = st.columns(5)

# Determine concrete checking status
stem_pass = M_actual_stem <= M_u_stem_allow if M_u_stem_allow > 0 else False
toe_pass = M_actual_toe <= M_u_ftg_sag_allow if M_u_ftg_sag_allow > 0 else False
heel_pass = M_actual_heel <= M_u_ftg_hog_allow if M_u_ftg_hog_allow > 0 else False
concrete_pass = stem_pass and toe_pass and heel_pass

with kpi_cols[0]:
    with st.container(border=True):
        if FS_ot_global is not None:
            status_ot = "PASS (Min 1.5)" if FS_ot_global >= 1.5 else "FAIL (Min 1.5)"
            delta_ot = "normal" if FS_ot_global >= 1.5 else "inverse"
            st.metric(label="Global Overturning FS", value=f"{FS_ot_global:.3f}", delta=status_ot, delta_color=delta_ot)
        else:
            st.metric(label="Global Overturning FS", value="N/A")

with kpi_cols[1]:
    with st.container(border=True):
        if FS_slide_global is not None:
            status_slide = "PASS (Min 1.5)" if FS_slide_global >= 1.5 else "FAIL (Min 1.5)"
            delta_slide = "normal" if FS_slide_global >= 1.5 else "inverse"
            st.metric(label="Global Sliding FS", value=f"{FS_slide_global:.3f}", delta=status_slide, delta_color=delta_slide)
        else:
            st.metric(label="Global Sliding FS", value="N/A")

with kpi_cols[2]:
    with st.container(border=True):
        if FS_ot_stem is not None:
            status_stem = "PASS (Min 1.5)" if FS_ot_stem >= 1.5 else "FAIL (Min 1.5)"
            delta_stem = "normal" if FS_ot_stem >= 1.5 else "inverse"
            st.metric(label="Stem Junction Overturning FS", value=f"{FS_ot_stem:.3f}", delta=status_stem, delta_color=delta_stem)
        else:
            st.metric(label="Stem Junction Overturning FS", value="N/A")

if enable_bearing:
    with kpi_cols[3]:
        with st.container(border=True):
            status_bearing = f"PASS (Min {FS_bearing:.2f})" if bearing_pass else f"FAIL (Min {FS_bearing:.2f})"
            delta_bearing = "normal" if bearing_pass else "inverse"
            bearing_fs_str = f"{FS_actual_bearing:.3f}" if FS_actual_bearing != float('inf') else "∞"
            st.metric(label="Bearing Capacity FS", value=bearing_fs_str, delta=status_bearing, delta_color=delta_bearing)
    disp_col_idx = 4
    concrete_col_idx = 5
else:
    disp_col_idx = 3
    concrete_col_idx = 4

with kpi_cols[disp_col_idx]:
    with st.container(border=True):
        st.metric(label="Max Nodal Disp (X / Y)", value=f"{max_disp_x*1000:.2f} / {max_disp_y*1000:.2f} mm", delta="FEM Peak", delta_color="off")

with kpi_cols[concrete_col_idx]:
    with st.container(border=True):
        st.metric(label="Concrete Strength", 
                  value="PASS" if concrete_pass else "FAIL", 
                  delta="Design Safe" if concrete_pass else "Overstressed",
                  delta_color="normal" if concrete_pass else "inverse")

# 3. Two-Column Dashboard Content Grid
col_left, col_right = st.columns([1, 1])

# --- LEFT COLUMN: MODEL & LOADS ---
with col_left:
    
    # Card A: Model Geometry & FE Mesh
    with st.container(border=True):
        st.subheader("📐 Model Geometry & FEM Mesh")
        st.markdown('<div class="section-desc">Dimensioned drawing of retaining wall components, soil blocks, water levels, and centroids.</div>', unsafe_allow_html=True)
        
        # Draw the model dynamically using drawsvg
        try:
            mult = 1000
            # Dynamically compute drawing height to prevent clipping when soil or water is high
            max_height = max(Hw, h_soil, h_soil_toe, Hwtr, Hwtr_front) + h_ftg
            d = draw.Drawing((ftg + 1.8)*mult, (max_height + 2.8)*mult, origin='bottom-left')
            
            # Draw soil - Dry soil (above water table)
            if h_soil > Hwtr:
                d.append(draw.Lines((toe+bot_wall)*mult, -(h_wet_height+h_ftg)*mult,
                                        ftg*mult, -(h_wet_height+h_ftg)*mult,
                                        ftg*mult, -(h_soil+h_ftg)*mult,
                                        (toe+bot_wall)*mult, -(h_soil+h_ftg)*mult,
                                        close=True,
                                fill='#5BC2A5',
                                stroke='black',
                                stroke_width=3))

            # Saturated soil (below water table)
            if h_wet_height > 0:
                d.append(draw.Lines(ftg*mult, -h_ftg*mult,
                                        ftg*mult, -(h_wet_height+h_ftg)*mult,
                                        (toe+bot_wall)*mult, -(h_wet_height+h_ftg)*mult,
                                        (toe+bot_wall)*mult, -h_ftg*mult,
                                        close=True,
                                fill='#A6F527',
                                stroke='black',
                                stroke_width=3))

            # Soil in front of toe
            if h_soil_toe > 0:
                d.append(draw.Lines(0, -h_ftg*mult,
                                        toe*mult, -h_ftg*mult,
                                        toe*mult, -(h_ftg + h_soil_toe)*mult,
                                        0, -(h_ftg + h_soil_toe)*mult,
                                        close=True,
                                fill='#5BC2A5',
                                stroke='black',
                                stroke_width=3))

            # Water level lines
            # Front water level
            if Hwtr_front > 0:
                d.append(draw.Line(0, -(h_ftg + Hwtr_front)*mult,
                                        toe*mult, -(h_ftg + Hwtr_front)*mult,
                                        stroke='blue', stroke_width=15, stroke_dasharray='50,50'))

            # Back water level
            if Hwtr > 0:
                d.append(draw.Line((toe+bot_wall)*mult, -(h_ftg + Hwtr)*mult,
                                        ftg*mult, -(h_ftg + Hwtr)*mult,
                                        stroke='blue', stroke_width=15, stroke_dasharray='50,50'))

            # Draw Retaining wall
            d.append(draw.Lines(0.0,  0.0,
                                    0.0,  -h_ftg*mult,
                                    toe*mult, -h_ftg*mult,
                                    (toe+(taper/4))*mult, -(h_ftg+(Hw/4))*mult,
                                    (toe+(taper/2))*mult, -(h_ftg+(Hw/2))*mult,
                                    (toe+(taper*(3/4)))*mult, -(h_ftg+(Hw*(3/4)))*mult,
                                    (toe+taper)*mult, -(h_ftg+Hw)*mult,
                                    (toe+bot_wall)*mult, -(h_ftg+Hw)*mult,
                                    (toe+bot_wall)*mult, -(h_ftg+(Hw*(3/4)))*mult,
                                    (toe+bot_wall)*mult, -(h_ftg+(Hw/2))*mult,
                                    (toe+bot_wall)*mult, -(h_ftg+(Hw/4))*mult,
                                    (toe+bot_wall)*mult, -h_ftg*mult,
                                    ftg*mult, -h_ftg*mult,
                                    ftg*mult, 0.0,
                                    close=True,
                            fill='#eeee00',
                            stroke='black',
                            stroke_width=50))

            # Vertical dimensions
            x_dim_v1 = (ftg + 0.3) * mult
            x_dim_v2 = (ftg + 0.9) * mult

            d.append(draw.Lines(ftg*mult, 0, x_dim_v2, 0, stroke='gray', stroke_width=10))
            d.append(draw.Lines(ftg*mult, -h_ftg*mult, x_dim_v1, -h_ftg*mult, stroke='gray', stroke_width=10))
            d.append(draw.Lines((toe+bot_wall)*mult, -(h_ftg+Hw)*mult, x_dim_v2, -(h_ftg+Hw)*mult, stroke='gray', stroke_width=10))

            # 1. Footing thickness (vertical)
            d.append(draw.Lines(x_dim_v1, 0, x_dim_v1, -h_ftg*mult, stroke='black', stroke_width=20))
            d.append(draw.Circle(x_dim_v1, 0, 0.06*mult))
            d.append(draw.Circle(x_dim_v1, -h_ftg*mult, 0.06*mult))
            d.append(draw.Text(f"{h_ftg} m", 0.2*mult, x_dim_v1 + (0.05*mult), -h_ftg*mult/2, fill='black'))

            # 2. Wall height (vertical)
            d.append(draw.Lines(x_dim_v1, -h_ftg*mult, x_dim_v1, -(h_ftg+Hw)*mult, stroke='black', stroke_width=20))
            d.append(draw.Circle(x_dim_v1, -h_ftg*mult, 0.06*mult))
            d.append(draw.Circle(x_dim_v1, -(h_ftg+Hw)*mult, 0.06*mult))
            d.append(draw.Text(f"{Hw} m", 0.2*mult, x_dim_v1 + (0.05*mult), -(h_ftg + Hw/2)*mult, fill='black'))

            # Total height (vertical)
            d.append(draw.Lines(x_dim_v2, 0, x_dim_v2, -(h_ftg+Hw)*mult, stroke='black', stroke_width=30))
            d.append(draw.Circle(x_dim_v2, 0, 0.08*mult))
            d.append(draw.Circle(x_dim_v2, -(h_ftg+Hw)*mult, 0.08*mult))
            d.append(draw.Text(f"{Hw+h_ftg} m", 0.25*mult, x_dim_v2 + (0.1*mult), -(h_ftg+Hw)*mult/2, fill='black'))

            # Horizontal dimensions
            y_dim_h1 = -(h_ftg + Hw + 0.3) * mult
            y_dim_h2 = -(h_ftg + Hw + 0.8) * mult
            y_dim_h3 = -(h_ftg + Hw + 1.4) * mult

            d.append(draw.Lines(0, 0, 0, y_dim_h3, stroke='gray', stroke_width=10))
            d.append(draw.Lines(toe*mult, -h_ftg*mult, toe*mult, y_dim_h2, stroke='gray', stroke_width=10))
            d.append(draw.Lines((toe+taper)*mult, -(h_ftg+Hw)*mult, (toe+taper)*mult, y_dim_h1, stroke='gray', stroke_width=10))
            d.append(draw.Lines((toe+bot_wall)*mult, -(h_ftg+Hw)*mult, (toe+bot_wall)*mult, y_dim_h2, stroke='gray', stroke_width=10))
            d.append(draw.Lines(ftg*mult, 0, ftg*mult, y_dim_h3, stroke='gray', stroke_width=10))

            # 3. Top wall thickness
            x_top_1 = (toe + taper) * mult
            x_top_2 = (toe + bot_wall) * mult
            d.append(draw.Lines(x_top_1, y_dim_h1, x_top_2, y_dim_h1, stroke='black', stroke_width=20))
            d.append(draw.Circle(x_top_1, y_dim_h1, 0.06*mult))
            d.append(draw.Circle(x_top_2, y_dim_h1, 0.06*mult))
            d.append(draw.Text(f"{top_wall} m", 0.2*mult, (x_top_1+x_top_2)/2, y_dim_h1 - (0.05*mult), fill='black', text_anchor='middle'))

            # 4. Toe length
            x_toe_1 = 0
            x_toe_2 = toe * mult
            d.append(draw.Lines(x_toe_1, y_dim_h2, x_toe_2, y_dim_h2, stroke='black', stroke_width=20))
            d.append(draw.Circle(x_toe_1, y_dim_h2, 0.06*mult))
            d.append(draw.Circle(x_toe_2, y_dim_h2, 0.06*mult))
            d.append(draw.Text(f"{toe} m", 0.2*mult, (x_toe_1+x_toe_2)/2, y_dim_h2 - (0.05*mult), fill='black', text_anchor='middle'))

            # 5. Bottom wall thickness
            x_bot_1 = toe * mult
            x_bot_2 = (toe + bot_wall) * mult
            d.append(draw.Lines(x_bot_1, y_dim_h2, x_bot_2, y_dim_h2, stroke='black', stroke_width=20))
            d.append(draw.Circle(x_bot_1, y_dim_h2, 0.06*mult))
            d.append(draw.Circle(x_bot_2, y_dim_h2, 0.06*mult))
            d.append(draw.Text(f"{bot_wall} m", 0.2*mult, (x_bot_1+x_bot_2)/2, y_dim_h2 - (0.05*mult), fill='black', text_anchor='middle'))

            # 6. Heel length
            x_heel_1 = (toe + bot_wall) * mult
            x_heel_2 = ftg * mult
            d.append(draw.Lines(x_heel_1, y_dim_h2, x_heel_2, y_dim_h2, stroke='black', stroke_width=20))
            d.append(draw.Circle(x_heel_1, y_dim_h2, 0.06*mult))
            d.append(draw.Circle(x_heel_2, y_dim_h2, 0.06*mult))
            d.append(draw.Text(f"{heel} m", 0.2*mult, (x_heel_1+x_heel_2)/2, y_dim_h2 - (0.05*mult), fill='black', text_anchor='middle'))

            # Total footing length
            d.append(draw.Lines(0, y_dim_h3, ftg*mult, y_dim_h3, stroke='black', stroke_width=30))
            d.append(draw.Circle(0, y_dim_h3, 0.08*mult))
            d.append(draw.Circle(ftg*mult, y_dim_h3, 0.08*mult))
            d.append(draw.Text(f"{ftg} m", 0.25*mult, (ftg*mult)/2, y_dim_h3 - (0.1*mult), fill='black', text_anchor='middle'))

            d.set_render_size(600, 400)

            # add centroid markers (red) for base, stem, and soil
            r = 0.03 * mult
            def add_x_lines(poly):
                if poly and len(poly) > 0:
                    tl = min(poly, key=lambda p: p[0] + p[1])
                    br = max(poly, key=lambda p: p[0] + p[1])
                    tr = max(poly, key=lambda p: p[0] - p[1])
                    bl = min(poly, key=lambda p: p[0] - p[1])
                    d.append(draw.Lines(tl[0] * mult, tl[1] * mult, br[0] * mult, br[1] * mult, stroke='black', stroke_dasharray='50,50', stroke_width=10))
                    d.append(draw.Lines(bl[0] * mult, bl[1] * mult, tr[0] * mult, tr[1] * mult, stroke='black', stroke_dasharray='50,50', stroke_width=10))

            if stem_centroid[0] is not None:
                sx, sy = stem_centroid
                add_x_lines(poly=stem_poly)
                d.append(draw.Circle(sx * mult, sy * mult, r * 1.4, fill='red', stroke='black', stroke_width=2))
            if base_centroid[0] is not None:
                bx, by = base_centroid
                add_x_lines(poly=base_poly)
                d.append(draw.Circle(bx * mult, by * mult, r * 1.4, fill='red', stroke='black', stroke_width=2))
            if 'soil_dry_centroid' in locals() and soil_dry_centroid[0] is not None:
                sd_x, sd_y = soil_dry_centroid
                add_x_lines(poly=soil_dry_poly)
                d.append(draw.Circle(sd_x * mult, sd_y * mult, r * 1.4, fill='red', stroke='black', stroke_width=2))
            if 'soil_wet_centroid' in locals() and soil_wet_centroid[0] is not None:
                sw_x, sw_y = soil_wet_centroid
                add_x_lines(poly=soil_wet_poly)
                d.append(draw.Circle(sw_x * mult, sw_y * mult, r * 1.4, fill='red', stroke='black', stroke_width=2))
                
            d.save_png('Rtwall.png')
        except Exception as e:
            st.error(f"Error drawing SVG model: {e}")

        # Display the SVG directly
        try:
            st.image(d.as_svg(), use_container_width=True)
        except Exception as e_svg:
            st.error(f"Error displaying SVG directly: {e_svg}")
            # Fallback to static png image
            st.image('Rtwall.png', use_container_width=True)
        
        # Expandable OpenSees Node Plot
        with st.expander("🔍 View OpenSees Finite Element Mesh Plot"):
            plt.close('all')
            opsv.plot_model("nodes", "elements")
            fig_model = plt.gcf()
            fig_model.set_size_inches(6, 4.5)
            if getattr(fig_model, 'axes', None):
                ax = fig_model.axes[0]
                try:
                    ax.relim()
                    ax.autoscale_view()
                except Exception:
                    pass
            fig_model.tight_layout()
            st.pyplot(fig_model)
            plt.close(fig_model)

    # Card A.2: SNI Geotechnical Geometry Verification
    with st.container(border=True):
        st.subheader("🛡️ SNI Geotechnical Geometry Verification Check")
        st.markdown('<div class="section-desc">Verification of typical dimensions based on Indonesian Geotechnical Standard (SNI Perencanaan Geoteknik - Gambar 35).</div>', unsafe_allow_html=True)

        try:
            # 1. Proportions checks
            H_tot = Hw + h_ftg
            top_wall_ok = top_wall >= 0.3
            bot_wall_ok = bot_wall >= 0.1 * H_tot
            slope_ok = (taper / Hw) >= (1.0 / 48.0) if Hw > 0 else True
            ftg_ok = (0.4 * H_tot) <= ftg <= (0.7 * H_tot)
            h_ftg_ok = (H_tot / 12.0) <= h_ftg <= (H_tot / 10.0)
            toe_ok = toe >= (ftg / 3.0)
            
            # 2. Draw verification model using drawsvg
            mult = 1000
            y_shift_sni = 1.5 * mult
            x_shift_sni = 0.8 * mult
            
            d_sni = draw.Drawing((ftg + 2.2)*mult, (H_tot + 1.8)*mult, origin='bottom-left')
            
            # Helper to draw a CAD-style dimension label directly on or near the line
            def draw_dim_cad_sni(d, x1, y1, x2, y2, text_label, value_str, req_str, is_ok, is_vertical=False, text_pos='right_or_top', text_y_override=None):
                color = '#27AE60' if is_ok else '#C0392B'
                status_str = "[PASS]" if is_ok else "[FAIL]"
                
                # 1. Main dimension line
                d.append(draw.Line(x1, y1, x2, y2, stroke=color, stroke_width=10))
                
                # 2. Tick lines at ends
                tick = 0.08 * mult
                if is_vertical:
                    d.append(draw.Line(x1 - tick, y1, x1 + tick, y1, stroke=color, stroke_width=10))
                    d.append(draw.Line(x2 - tick, y2, x2 + tick, y2, stroke=color, stroke_width=10))
                    
                    # Aligned text label to the right
                    tx = x1 + 0.15 * mult
                    ty = (y1 + y2) / 2
                    d.append(draw.Text(f"{text_label} = {value_str}", 0.18*mult, tx, ty + 0.1*mult, fill=color, font_weight='bold'))
                    d.append(draw.Text(f"Req: {req_str} {status_str}", 0.14*mult, tx, ty - 0.1*mult, fill=color))
                else:
                    d.append(draw.Line(x1, y1 - tick, x1, y1 + tick, stroke=color, stroke_width=10))
                    d.append(draw.Line(x2, y2 - tick, x2, y2 + tick, stroke=color, stroke_width=10))
                    
                    # Aligned text label above or below
                    tx = (x1 + x2) / 2
                    if text_y_override is not None:
                        ty = text_y_override
                        d.append(draw.Text(f"{text_label} = {value_str}", 0.18*mult, tx, ty, fill=color, font_weight='bold', text_anchor='middle'))
                        d.append(draw.Text(f"Req: {req_str} {status_str}", 0.14*mult, tx, ty + 0.20*mult, fill=color, text_anchor='middle'))
                    else:
                        if text_pos == 'top':
                            ty = y1 + 0.1 * mult
                            d.append(draw.Text(f"{text_label} = {value_str}", 0.18*mult, tx, ty + 0.16*mult, fill=color, font_weight='bold', text_anchor='middle'))
                            d.append(draw.Text(f"Req: {req_str} {status_str}", 0.14*mult, tx, ty, fill=color, text_anchor='middle'))
                        else:
                            ty = y1 - 0.18 * mult
                            d.append(draw.Text(f"{text_label} = {value_str}", 0.18*mult, tx, ty, fill=color, font_weight='bold', text_anchor='middle'))
                            d.append(draw.Text(f"Req: {req_str} {status_str}", 0.14*mult, tx, ty - 0.15*mult, fill=color, text_anchor='middle'))

            # Draw retaining wall structure (yellow)
            d_sni.append(draw.Lines(x_shift_sni,  -y_shift_sni,
                                    x_shift_sni,  -h_ftg*mult - y_shift_sni,
                                    toe*mult + x_shift_sni, -h_ftg*mult - y_shift_sni,
                                    (toe+(taper/4))*mult + x_shift_sni, -(h_ftg+(Hw/4))*mult - y_shift_sni,
                                    (toe+(taper/2))*mult + x_shift_sni, -(h_ftg+(Hw/2))*mult - y_shift_sni,
                                    (toe+(taper*(3/4)))*mult + x_shift_sni, -(h_ftg+(Hw*(3/4)))*mult - y_shift_sni,
                                    (toe+taper)*mult + x_shift_sni, -(h_ftg+Hw)*mult - y_shift_sni,
                                    (toe+bot_wall)*mult + x_shift_sni, -(h_ftg+Hw)*mult - y_shift_sni,
                                    (toe+bot_wall)*mult + x_shift_sni, -(h_ftg+(Hw*(3/4)))*mult - y_shift_sni,
                                    (toe+bot_wall)*mult + x_shift_sni, -(h_ftg+(Hw/2))*mult - y_shift_sni,
                                    (toe+bot_wall)*mult + x_shift_sni, -(h_ftg+(Hw/4))*mult - y_shift_sni,
                                    (toe+bot_wall)*mult + x_shift_sni, -h_ftg*mult - y_shift_sni,
                                    ftg*mult + x_shift_sni, -h_ftg*mult - y_shift_sni,
                                    ftg*mult + x_shift_sni, -y_shift_sni,
                                    close=True,
                                    fill='#eeee00',
                                    stroke='black',
                                    stroke_width=25))

            # 1. Total Height H (vertical on right)
            draw_dim_cad_sni(d_sni, (ftg + 0.15)*mult + x_shift_sni, -y_shift_sni,
                                     (ftg + 0.15)*mult + x_shift_sni, -H_tot*mult - y_shift_sni,
                                     "Height H", f"{H_tot:.2f} m", "Reference", True, is_vertical=True)

            # 2. Footing Thickness D (vertical on right)
            draw_dim_cad_sni(d_sni, (ftg + 0.85)*mult + x_shift_sni, -y_shift_sni,
                                     (ftg + 0.85)*mult + x_shift_sni, -h_ftg*mult - y_shift_sni,
                                     "Thickness D", f"{h_ftg:.2f} m", f"{H_tot/12:.2f}~{H_tot/10:.2f} m", h_ftg_ok, is_vertical=True)

            # Aligned Y coordinate for toe length, Width B, and base wall dimension lines (below the bottom of footing -1500)
            y_dim_bottom = -y_shift_sni + 0.15*mult

            # 3. Footing Width B (at bottom) - Text level 3 (Width B text below base wall text)
            draw_dim_cad_sni(d_sni, x_shift_sni, y_dim_bottom,
                                     ftg*mult + x_shift_sni, y_dim_bottom,
                                     "Width B", f"{ftg:.2f} m", f"{0.4*H_tot:.2f}~{0.7*H_tot:.2f} m", ftg_ok, is_vertical=False, text_y_override=y_dim_bottom + 1.0*mult)

            # 4. Top Wall thickness (above stem top)
            draw_dim_cad_sni(d_sni, (toe + taper)*mult + x_shift_sni, -(h_ftg + Hw + 0.15)*mult - y_shift_sni,
                                     (toe + bot_wall)*mult + x_shift_sni, -(h_ftg + Hw + 0.15)*mult - y_shift_sni,
                                     "Top Wall", f"{top_wall:.2f} m", ">= 0.30 m", top_wall_ok, is_vertical=False, text_pos='top')

            # 5. Bottom Wall thickness (Base Wall) - Text level 2 (base wall text below toe length text)
            draw_dim_cad_sni(d_sni, toe*mult + x_shift_sni, y_dim_bottom,
                                     (toe + bot_wall)*mult + x_shift_sni, y_dim_bottom,
                                     "Base Wall", f"{bot_wall:.2f} m", f">={0.1*H_tot:.2f} m", bot_wall_ok, is_vertical=False, text_y_override=y_dim_bottom + 0.5*mult)

            # 6. Toe length - Text level 1 (first toe length text closest to dimension line)
            draw_dim_cad_sni(d_sni, x_shift_sni, y_dim_bottom,
                                     toe*mult + x_shift_sni, y_dim_bottom,
                                     "Toe Length", f"{toe:.2f} m", f">={ftg/3:.2f} m", toe_ok, is_vertical=False, text_y_override=y_dim_bottom + 0.0*mult)

            # 7. Front Slope Batter
            slope_val = taper / Hw if Hw > 0 else 0
            slope_str = f"1 : {1.0/slope_val:.1f}" if slope_val > 0 else "Vertical"
            draw_dim_cad_sni(d_sni, (toe + taper/2)*mult + x_shift_sni - 0.15*mult, -(h_ftg + Hw/2)*mult - y_shift_sni,
                                     (toe + taper/2)*mult + x_shift_sni + 0.15*mult, -(h_ftg + Hw/2)*mult - y_shift_sni,
                                     "Batter", slope_str, ">= 1:48", slope_ok, is_vertical=False, text_pos='top')

            d_sni.set_render_size(800, 550)
            st.image(d_sni.as_svg(), use_container_width=True)

            # Proportions checklist table below
            df_checks = pd.DataFrame([
                {
                    "Parameter Name": "Height (H)",
                    "Parameter Symbol": "H",
                    "Input Value": f"{H_tot:.2f} m",
                    "SNI Requirement": "Reference base",
                    "Status": "INFO"
                },
                {
                    "Parameter Name": "Top Wall Thickness",
                    "Parameter Symbol": "top_wall",
                    "Input Value": f"{top_wall:.2f} m",
                    "SNI Requirement": "≥ 0.30 m",
                    "Status": "PASS" if top_wall_ok else "FAIL"
                },
                {
                    "Parameter Name": "Base Stem Thickness",
                    "Parameter Symbol": "bot_wall",
                    "Input Value": f"{bot_wall:.2f} m",
                    "SNI Requirement": f"≥ 0.1 H ({0.1*H_tot:.2f} m)",
                    "Status": "PASS" if bot_wall_ok else "FAIL"
                },
                {
                    "Parameter Name": "Front Face Batter Slope",
                    "Parameter Symbol": "slope",
                    "Input Value": f"{slope_str} ({slope_val:.4f})",
                    "SNI Requirement": "≥ 1:48 (0.0208)",
                    "Status": "PASS" if slope_ok else "FAIL"
                },
                {
                    "Parameter Name": "Footing Width (B)",
                    "Parameter Symbol": "ftg",
                    "Input Value": f"{ftg:.2f} m",
                    "SNI Requirement": f"0.4 H ~ 0.7 H ({0.4*H_tot:.2f} ~ {0.7*H_tot:.2f} m)",
                    "Status": "PASS" if ftg_ok else "FAIL"
                },
                {
                    "Parameter Name": "Footing Thickness (D)",
                    "Parameter Symbol": "h_ftg",
                    "Input Value": f"{h_ftg:.2f} m",
                    "SNI Requirement": f"H/12 ~ H/10 ({H_tot/12:.2f} ~ {H_tot/10:.2f} m)",
                    "Status": "PASS" if h_ftg_ok else "FAIL"
                },
                {
                    "Parameter Name": "Toe Slab Length",
                    "Parameter Symbol": "toe",
                    "Input Value": f"{toe:.2f} m",
                    "SNI Requirement": f"≥ B/3 ({ftg/3:.2f} m)",
                    "Status": "PASS" if toe_ok else "FAIL"
                }
            ])

            # Style helper
            def style_status(val):
                if val == "PASS":
                    return "color: #27AE60; font-weight: bold;"
                elif val == "FAIL":
                    return "color: #C0392B; font-weight: bold;"
                return "color: #2980B9; font-weight: bold;"

            try:
                styled_df = df_checks.style.map(style_status, subset=["Status"])
            except AttributeError:
                styled_df = df_checks.style.applymap(style_status, subset=["Status"])
            st.table(styled_df)

            # Combined status alert
            all_sni_ok = all([top_wall_ok, bot_wall_ok, slope_ok, ftg_ok, h_ftg_ok, toe_ok])
            if all_sni_ok:
                st.success("🎉 **Success:** All wall geometry parameters meet the typical proportions recommended by **SNI Perencanaan Geoteknik** (Gambar 35b).")
            else:
                st.warning("⚠️ **Warning:** One or more geometry parameters do not satisfy the typical proportions of **SNI Perencanaan Geoteknik** (Gambar 35b). Please adjust the input dimensions to comply.")

        except Exception as ex:
            st.error(f"Error executing SNI Geotechnical verification: {ex}")

    # Card B: Applied Loads & Pressures
    with st.container(border=True):
        st.subheader("🌊 Applied Loads & Pressure Distributions")
        st.markdown('<div class="section-desc">Visual and tabular representation of boundary loading, earth/water pressures, and surcharges.</div>', unsafe_allow_html=True)
        
        # Dropdown list to change type of load displayed
        selected_load_plot = st.selectbox(
            "Select Load Type to Plot",
            [
                "total",
                "soil lateral",
                "soil vertical",
                "hydrostatic",
                "uplift",
                "seismic"
            ],
            index=0,
            key="selected_load_plot"
        )

        # Draw the model dynamically using drawsvg for loads
        try:
            mult = 1000
            y_shift = 0.8 * mult
            x_shift = 0.8 * mult
            max_height = max(Hw, h_soil, h_soil_toe, Hwtr, Hwtr_front) + h_ftg
            
            # Initialize drawing
            d_load = draw.Drawing((ftg + 3.0)*mult, (max_height + 3.0)*mult, origin='bottom-left')
            
            # Helper for arrows
            def draw_arrow_head(d, x1, y1, x2, y2, color='red', stroke_width=15, head_len=80, head_width=50):
                d.append(draw.Line(x1, y1, x2, y2, stroke=color, stroke_width=stroke_width))
                # Calculate arrowhead
                dx = x2 - x1
                dy = y2 - y1
                L = np.hypot(dx, dy)
                if L < 1e-6:
                    return
                ux = dx / L
                uy = dy / L
                # Perp vector
                px = -uy
                py = ux
                # Arrow head points
                ax = x2 - head_len * ux
                ay = y2 - head_len * uy
                
                p1x = ax + head_width * px
                p1y = ay + head_width * py
                p2x = ax - head_width * px
                p2y = ay - head_width * py
                
                d.append(draw.Lines(x2, y2, p1x, p1y, p2x, p2y, close=True, fill=color, stroke=color))

            # 1. Background Soil outlines/shapes (very light)
            # Soil dry
            if h_soil > Hwtr:
                d_load.append(draw.Lines((toe+bot_wall)*mult + x_shift, -(h_wet_height+h_ftg)*mult - y_shift,
                                     ftg*mult + x_shift, -(h_wet_height+h_ftg)*mult - y_shift,
                                     ftg*mult + x_shift, -(h_soil+h_ftg)*mult - y_shift,
                                     (toe+bot_wall)*mult + x_shift, -(h_soil+h_ftg)*mult - y_shift,
                                     close=True, fill='#5BC2A5', fill_opacity=0.06, stroke='#bdc3c7', stroke_width=2, stroke_dasharray='10,10'))
            # Soil wet
            if h_wet_height > 0:
                d_load.append(draw.Lines(ftg*mult + x_shift, -h_ftg*mult - y_shift,
                                     ftg*mult + x_shift, -(h_wet_height+h_ftg)*mult - y_shift,
                                     (toe+bot_wall)*mult + x_shift, -(h_wet_height+h_ftg)*mult - y_shift,
                                     (toe+bot_wall)*mult + x_shift, -h_ftg*mult - y_shift,
                                     close=True, fill='#A6F527', fill_opacity=0.06, stroke='#bdc3c7', stroke_width=2, stroke_dasharray='10,10'))
            # Soil toe
            if h_soil_toe > 0:
                d_load.append(draw.Lines(x_shift, -h_ftg*mult - y_shift,
                                     toe*mult + x_shift, -h_ftg*mult - y_shift,
                                     toe*mult + x_shift, -(h_ftg + h_soil_toe)*mult - y_shift,
                                     x_shift, -(h_ftg + h_soil_toe)*mult - y_shift,
                                     close=True, fill='#5BC2A5', fill_opacity=0.06, stroke='#bdc3c7', stroke_width=2, stroke_dasharray='10,10'))

            # Water level lines (dashed)
            if Hwtr_front > 0:
                d_load.append(draw.Line(x_shift, -(h_ftg + Hwtr_front)*mult - y_shift,
                                        toe*mult + x_shift, -(h_ftg + Hwtr_front)*mult - y_shift,
                                        stroke='blue', stroke_width=10, stroke_dasharray='30,30'))
            if Hwtr > 0:
                d_load.append(draw.Line((toe+bot_wall)*mult + x_shift, -(h_ftg + Hwtr)*mult - y_shift,
                                        ftg*mult + x_shift, -(h_ftg + Hwtr)*mult - y_shift,
                                        stroke='blue', stroke_width=10, stroke_dasharray='30,30'))

            # 2. Draw yellow wall
            d_load.append(draw.Lines(x_shift,  -y_shift,
                                x_shift,  -h_ftg*mult - y_shift,
                                toe*mult + x_shift, -h_ftg*mult - y_shift,
                                (toe+(taper/4))*mult + x_shift, -(h_ftg+(Hw/4))*mult - y_shift,
                                (toe+(taper/2))*mult + x_shift, -(h_ftg+(Hw/2))*mult - y_shift,
                                (toe+(taper*(3/4)))*mult + x_shift, -(h_ftg+(Hw*(3/4)))*mult - y_shift,
                                (toe+taper)*mult + x_shift, -(h_ftg+Hw)*mult - y_shift,
                                (toe+bot_wall)*mult + x_shift, -(h_ftg+Hw)*mult - y_shift,
                                (toe+bot_wall)*mult + x_shift, -(h_ftg+(Hw*(3/4)))*mult - y_shift,
                                (toe+bot_wall)*mult + x_shift, -(h_ftg+(Hw/2))*mult - y_shift,
                                (toe+bot_wall)*mult + x_shift, -(h_ftg+(Hw/4))*mult - y_shift,
                                (toe+bot_wall)*mult + x_shift, -h_ftg*mult - y_shift,
                                ftg*mult + x_shift, -h_ftg*mult - y_shift,
                                ftg*mult + x_shift, -y_shift,
                                close=True,
                                fill='#eeee00',
                                stroke='black',
                                stroke_width=30))

            # 3. Dynamic loads
            if selected_load_plot in ['soil lateral', 'total']:
                # Soil lateral pressure profile: active soil + surcharge
                if surcharge_type == 'Strip Load':
                    q_top = stresses_stripload_retainingwall_local(q, width_surcharge, offset_surcharge, h_soil, max(1e-5, h_soil - Hw))['delta sigma x [kPa]']
                    q_bot = stresses_stripload_retainingwall_local(q, width_surcharge, offset_surcharge, h_soil, h_soil)['delta sigma x [kPa]']
                else:
                    q_top = Ka * q
                    q_bot = Ka * q
                p_top = Ka * (gamma_soil_dry * max(0.0, h_soil - Hw)) + q_top
                p_bot = Ka * (gamma_soil_dry * max(0.0, h_soil - Hwtr) + gamma_soil_wet * min(h_soil, Hwtr)) + q_bot
                
                # Scale: max pressure maps to 1.2 meters screen-wise
                p_max_lat = max(p_top, p_bot, 1.0)
                scale_lat = 1.2 * mult / p_max_lat
                
                # Draw shaded pressure block
                d_load.append(draw.Lines((toe+bot_wall)*mult + x_shift, -(h_ftg+Hw)*mult - y_shift,
                                    (toe+bot_wall)*mult + x_shift + p_top*scale_lat, -(h_ftg+Hw)*mult - y_shift,
                                    (toe+bot_wall)*mult + x_shift + p_bot*scale_lat, -h_ftg*mult - y_shift,
                                    (toe+bot_wall)*mult + x_shift, -h_ftg*mult - y_shift,
                                    close=True, fill='#E67E22', fill_opacity=0.3, stroke='#D35400', stroke_width=12))
                
                # Draw arrows (spaced evenly)
                n_arr = 5
                for i in range(n_arr):
                    frac = i / (n_arr - 1)
                    z_val = frac * Hw
                    hd = max(0.0, h_soil - max(z_val, Hwtr))
                    hw = max(0.0, min(h_soil, Hwtr) - z_val)
                    if surcharge_type == 'Strip Load':
                        q_z = stresses_stripload_retainingwall_local(q, width_surcharge, offset_surcharge, h_soil, max(1e-5, h_soil - z_val))['delta sigma x [kPa]']
                    else:
                        q_z = Ka * q
                    pz = Ka * (gamma_soil_dry * hd + gamma_soil_wet * hw) + q_z
                    y_z = -(h_ftg + z_val)*mult - y_shift
                    x_start = (toe + bot_wall)*mult + x_shift + pz * scale_lat
                    x_end = (toe + bot_wall)*mult + x_shift
                    if pz > 0:
                        draw_arrow_head(d_load, x_start, y_z, x_end, y_z, color='#D35400', stroke_width=10, head_len=60, head_width=35)
                
                # Label pressures
                d_load.append(draw.Text(f"{p_top:.2f} kPa", 0.20*mult, (toe+bot_wall)*mult + x_shift + p_top*scale_lat + 0.20*mult, -(h_ftg+Hw)*mult - y_shift, fill='black'))
                d_load.append(draw.Text(f"{p_bot:.2f} kPa", 0.20*mult, (toe+bot_wall)*mult + x_shift + p_bot*scale_lat + 0.20*mult, -h_ftg*mult - y_shift, fill='black'))

            if selected_load_plot in ['soil vertical', 'total']:
                # Soil vertical load on heel: weight + surcharge
                p_v = q_soil + q
                # Scale: max vertical pressure maps to 0.8 meters
                scale_v = 0.8 * mult / max(p_v, 1.0)
                
                # Draw shaded block above footing heel
                d_load.append(draw.Lines((toe+bot_wall)*mult + x_shift, -h_ftg*mult - y_shift,
                                    (toe+bot_wall)*mult + x_shift, -h_ftg*mult - y_shift - p_v*scale_v,
                                    ftg*mult + x_shift, -h_ftg*mult - y_shift - p_v*scale_v,
                                    ftg*mult + x_shift, -h_ftg*mult - y_shift,
                                    close=True, fill='#F39C12', fill_opacity=0.3, stroke='#E67E22', stroke_width=12))
                
                # Draw downward arrows
                n_arr = 5
                for i in range(n_arr):
                    frac = i / (n_arr - 1)
                    x_pos = (toe + bot_wall + frac * heel)*mult + x_shift
                    y_start = -h_ftg*mult - y_shift - p_v*scale_v
                    y_end = -h_ftg*mult - y_shift
                    draw_arrow_head(d_load, x_pos, y_start, x_pos, y_end, color='#E67E22', stroke_width=10, head_len=60, head_width=35)
                
                # Label pressure
                d_load.append(draw.Text(f"{p_v:.2f} kPa", 0.20*mult, (toe+bot_wall + heel/2)*mult + x_shift, -h_ftg*mult - y_shift - p_v*scale_v - 0.22*mult, fill='black', text_anchor='middle'))

            if selected_load_plot in ['hydrostatic', 'total']:
                # Hydrostatic water pressure
                p_w_back = gamma_w * Hwtr
                p_w_front = gamma_w * Hwtr_front
                p_max_w = max(p_w_back, p_w_front, 1.0)
                scale_w = 1.2 * mult / p_max_w
                
                # Back water pressure triangle
                if Hwtr > 0:
                    d_load.append(draw.Lines((toe+bot_wall)*mult + x_shift, -(h_ftg+Hwtr)*mult - y_shift,
                                        (toe+bot_wall)*mult + x_shift + p_w_back*scale_w, -h_ftg*mult - y_shift,
                                        (toe+bot_wall)*mult + x_shift, -h_ftg*mult - y_shift,
                                        close=True, fill='#3498DB', fill_opacity=0.3, stroke='#2980B9', stroke_width=12))
                    # Back water arrows
                    n_arr = 3
                    for i in range(n_arr):
                        frac = i / (n_arr - 1)
                        z_val = frac * Hwtr
                        pz_w = gamma_w * (Hwtr - z_val)
                        y_z = -(h_ftg + z_val)*mult - y_shift
                        x_start = (toe + bot_wall)*mult + x_shift + pz_w * scale_w
                        x_end = (toe + bot_wall)*mult + x_shift
                        if pz_w > 0:
                            draw_arrow_head(d_load, x_start, y_z, x_end, y_z, color='#2980B9', stroke_width=10, head_len=50, head_width=30)
                    
                    d_load.append(draw.Text(f"{p_w_back:.2f} kPa", 0.20*mult, (toe+bot_wall)*mult + x_shift + p_w_back*scale_w + 0.20*mult, -h_ftg*mult - y_shift, fill='black'))

                # Front water pressure triangle
                if Hwtr_front > 0:
                    d_load.append(draw.Lines(toe*mult + x_shift, -(h_ftg+Hwtr_front)*mult - y_shift,
                                        toe*mult + x_shift - p_w_front*scale_w, -h_ftg*mult - y_shift,
                                        toe*mult + x_shift, -h_ftg*mult - y_shift,
                                        close=True, fill='#3498DB', fill_opacity=0.3, stroke='#2980B9', stroke_width=12))
                    # Front water arrows (pointing right)
                    n_arr = 3
                    for i in range(n_arr):
                        frac = i / (n_arr - 1)
                        z_val = frac * Hwtr_front
                        pz_w = gamma_w * (Hwtr_front - z_val)
                        y_z = -(h_ftg + z_val)*mult - y_shift
                        x_start = toe*mult + x_shift - pz_w * scale_w
                        x_end = toe*mult + x_shift
                        if pz_w > 0:
                            draw_arrow_head(d_load, x_start, y_z, x_end, y_z, color='#2980B9', stroke_width=10, head_len=50, head_width=30)
                    
                    d_load.append(draw.Text(f"{p_w_front:.2f} kPa", 0.20*mult, toe*mult + x_shift - p_w_front*scale_w - 0.85*mult, -h_ftg*mult - y_shift, fill='black'))

            if selected_load_plot in ['uplift', 'total']:
                # Uplift pressure profile acting upwards under base
                p_upl = gamma_w * Hwtr
                if p_upl > 0:
                    # Scale: uplift pressure maps to 0.6 meters screen-wise
                    scale_u = 0.6 * mult / p_upl
                    
                    # Draw shaded uplift block below the base
                    d_load.append(draw.Lines(x_shift, -y_shift,
                                        x_shift, -y_shift + p_upl*scale_u,
                                        ftg*mult + x_shift, -y_shift + p_upl*scale_u,
                                        ftg*mult + x_shift, -y_shift,
                                        close=True, fill='#9B59B6', fill_opacity=0.3, stroke='#8E44AD', stroke_width=12))
                    
                    # Draw upward arrows
                    n_arr = 6
                    for i in range(n_arr):
                        frac = i / (n_arr - 1)
                        x_pos = frac * ftg * mult + x_shift
                        y_start = -y_shift + p_upl*scale_u
                        y_end = -y_shift
                        draw_arrow_head(d_load, x_pos, y_start, x_pos, y_end, color='#8E44AD', stroke_width=10, head_len=50, head_width=30)
                        
                    d_load.append(draw.Text(f"{p_upl:.2f} kPa", 0.20*mult, (ftg/2)*mult + x_shift, -y_shift + p_upl*scale_u + 0.25*mult, fill='black', text_anchor='middle'))

            if selected_load_plot in ['seismic', 'total']:
                # Pseudo-static horizontal seismic forces acting left at centroids
                if kh > 0:
                    # Draw arrows at centroids of stem, base, and backfill
                    # Stem centroid
                    if stem_centroid[0] is not None:
                        sx, sy = stem_centroid[0]*mult + x_shift, stem_centroid[1]*mult - y_shift
                        draw_arrow_head(d_load, sx + 0.8*mult, sy, sx, sy, color='red', stroke_width=15, head_len=80, head_width=45)
                        d_load.append(draw.Text("Seismic (stem)", 0.18*mult, sx + 0.95*mult, sy + 0.05*mult, fill='red', font_weight='bold'))
                    # Base centroid
                    if base_centroid[0] is not None:
                        bx, by = base_centroid[0]*mult + x_shift, base_centroid[1]*mult - y_shift
                        draw_arrow_head(d_load, bx + 0.8*mult, by, bx, by, color='red', stroke_width=15, head_len=80, head_width=45)
                        d_load.append(draw.Text("Seismic (base)", 0.18*mult, bx + 0.95*mult, by + 0.05*mult, fill='red', font_weight='bold'))
                    # Backfill centroids (dry / wet)
                    if 'soil_dry_centroid' in locals() and soil_dry_centroid[0] is not None:
                        sdx, sdy = soil_dry_centroid[0]*mult + x_shift, soil_dry_centroid[1]*mult - y_shift
                        draw_arrow_head(d_load, sdx + 0.8*mult, sdy, sdx, sdy, color='red', stroke_width=15, head_len=80, head_width=45)
                        d_load.append(draw.Text("Seismic (soil dry)", 0.18*mult, sdx + 0.95*mult, sdy + 0.05*mult, fill='red', font_weight='bold'))
                    if 'soil_wet_centroid' in locals() and soil_wet_centroid[0] is not None:
                        swx, swy = soil_wet_centroid[0]*mult + x_shift, soil_wet_centroid[1]*mult - y_shift
                        draw_arrow_head(d_load, swx + 0.8*mult, swy, swx, swy, color='red', stroke_width=15, head_len=80, head_width=45)
                        d_load.append(draw.Text("Seismic (soil wet)", 0.18*mult, swx + 0.95*mult, swy + 0.05*mult, fill='red', font_weight='bold'))

            if selected_load_plot == 'total':
                # Concrete weights downward at centroids
                if base_centroid[0] is not None:
                    bx, by = base_centroid[0]*mult + x_shift, base_centroid[1]*mult - y_shift
                    draw_arrow_head(d_load, bx, by - 0.8*mult, bx, by, color='blue', stroke_width=15, head_len=80, head_width=45)
                    d_load.append(draw.Text(f"W_base = {W_base:.1f} kN", 0.18*mult, bx, by - 1.0*mult, fill='blue', text_anchor='middle'))
                if stem_centroid[0] is not None:
                    sx, sy = stem_centroid[0]*mult + x_shift, stem_centroid[1]*mult - y_shift
                    draw_arrow_head(d_load, sx, sy - 0.8*mult, sx, sy, color='blue', stroke_width=15, head_len=80, head_width=45)
                    d_load.append(draw.Text(f"W_stem = {W_stem_rect:.1f} kN", 0.18*mult, sx, sy - 1.0*mult, fill='blue', text_anchor='middle'))

            d_load.set_render_size(800, 550)
            st.image(d_load.as_svg(), use_container_width=True)
            
        except Exception as e:
            st.error(f"Error drawing SVG load diagram: {e}")

        # Max Load Metrics Display under the diagram
        st.markdown("### 📊 Peak Loading Parameters")
        
        try:
            soil_lat_force = abs(sum(Fx_soil_q.values()))
        except Exception:
            soil_lat_force = 0.0
            
        try:
            water_lat_force = abs(sum(Fx_water_lat.values()))
        except Exception:
            water_lat_force = 0.0
            
        if selected_load_plot == 'soil lateral':
            if surcharge_type == 'Strip Load':
                q_top = stresses_stripload_retainingwall_local(q, width_surcharge, offset_surcharge, h_soil, max(1e-5, h_soil - Hw))['delta sigma x [kPa]']
                q_bot = stresses_stripload_retainingwall_local(q, width_surcharge, offset_surcharge, h_soil, h_soil)['delta sigma x [kPa]']
            else:
                q_top = Ka * q
                q_bot = Ka * q
            p_top = Ka * (gamma_soil_dry * max(0.0, h_soil - Hw)) + q_top
            p_bot = Ka * (gamma_soil_dry * max(0.0, h_soil - Hwtr) + gamma_soil_wet * min(h_soil, Hwtr)) + q_bot
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Max Lateral Pressure (Base)", f"{p_bot:.2f} kPa")
            with col2:
                st.metric("Min Lateral Pressure (Top)", f"{p_top:.2f} kPa")
            with col3:
                st.metric("Total Lateral Soil Force", f"{soil_lat_force:.2f} kN")
                
        elif selected_load_plot == 'soil vertical':
            if surcharge_type == 'Strip Load':
                n_pts = 10
                forces_pts = []
                for k in range(n_pts):
                    frac = (k + 0.5) / n_pts
                    d = frac * heel
                    x_rel = d - offset_surcharge
                    res = stresses_stripload_local(z=h_soil, x=x_rel, width=width_surcharge, imposedstress=q)
                    forces_pts.append(max(0.0, res['delta sigma z [kPa]']))
                avg_q_v = np.mean(forces_pts)
            else:
                avg_q_v = q
            p_v = q_soil + avg_q_v
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Avg Vertical Pressure on Heel", f"{p_v:.2f} kPa")
            with col2:
                st.metric("Total Soil Weight on Heel", f"{abs(soil_total):.2f} kN")
                
        elif selected_load_plot == 'hydrostatic':
            p_w_back = gamma_w * Hwtr
            p_w_front = gamma_w * Hwtr_front
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Max Water Pressure (Back)", f"{p_w_back:.2f} kPa")
            with col2:
                st.metric("Max Water Pressure (Front)", f"{p_w_front:.2f} kPa" if Hwtr_front > 0 else "0.00 kPa")
            with col3:
                st.metric("Total Hydrostatic Lateral Force", f"{water_lat_force:.2f} kN")
                
        elif selected_load_plot == 'uplift':
            p_uplift = gamma_w * Hwtr
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Uplift Pressure at Base", f"{p_uplift:.2f} kPa")
            with col2:
                st.metric("Total Base Uplift Force", f"{uplift_total:.2f} kN")
                
        elif selected_load_plot == 'seismic':
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Seismic Acceleration (kh)", f"{kh:.4f}")
            with col2:
                st.metric("Total Seismic Inertia Force", f"{F_seismic_pseudo:.2f} kN")
                
        elif selected_load_plot == 'total':
            col1, col2, col3 = st.columns(3)
            with col1:
                if surcharge_type == 'Strip Load':
                    q_bot = stresses_stripload_retainingwall_local(q, width_surcharge, offset_surcharge, h_soil, h_soil)['delta sigma x [kPa]']
                else:
                    q_bot = Ka * q
                p_bot = Ka * (gamma_soil_dry * max(0.0, h_soil - Hwtr) + gamma_soil_wet * min(h_soil, Hwtr)) + q_bot
                st.metric("Max Soil Lateral Pressure", f"{p_bot:.2f} kPa")
                st.metric("Soil Lateral Force", f"{soil_lat_force:.2f} kN")
            with col2:
                if surcharge_type == 'Strip Load':
                    n_pts = 10
                    forces_pts = []
                    for k in range(n_pts):
                        frac = (k + 0.5) / n_pts
                        d = frac * heel
                        x_rel = d - offset_surcharge
                        res = stresses_stripload_local(z=h_soil, x=x_rel, width=width_surcharge, imposedstress=q)
                        forces_pts.append(max(0.0, res['delta sigma z [kPa]']))
                    avg_q_v = np.mean(forces_pts)
                else:
                    avg_q_v = q
                p_v = q_soil + avg_q_v
                st.metric("Avg Soil Vertical Pressure", f"{p_v:.2f} kPa")
                st.metric("Soil Vertical Force", f"{abs(soil_total):.2f} kN")
            with col3:
                p_w_back = gamma_w * Hwtr
                st.metric("Max Hydrostatic Pressure", f"{p_w_back:.2f} kPa")
                st.metric("Hydrostatic Lateral Force", f"{water_lat_force:.2f} kN")
            
            st.markdown("---")
            col4, col5, col6 = st.columns(3)
            with col4:
                p_uplift = gamma_w * Hwtr
                st.metric("Uplift Pressure at Base", f"{p_uplift:.2f} kPa")
                st.metric("Total Uplift Force", f"{uplift_total:.2f} kN")
            with col5:
                st.metric("Seismic Accel (kh)", f"{kh:.4f}")
                st.metric("Seismic Inertial Force", f"{F_seismic_pseudo:.2f} kN")
            with col6:
                st.metric("Concrete Self-Weight", f"{W_conc:.2f} kN")
                st.metric("Net Downward Force (heel - uplift)", f"{abs(soil_total) + W_conc - uplift_total:.2f} kN")
            
        # Summary tables in expanders
        with st.expander("📊 Load Summary Tables (Node Forces)"):
            if not df_lateral.empty:
                st.markdown("**Summary Table — Lateral Loads (kN)**")
                st.table(df_lateral)
            if not df_vertical.empty:
                st.markdown("**Summary Table — Vertical Loads (kN)**")
                st.table(df_vertical)
                
        # Detailed Load calculations in expanders
        with st.expander("📝 Show Analytical Load Calculations"):
            st.markdown("#### Soil Lateral & Vertical Load")
            st.latex(r"K_a = \tan^2\left(45^\circ - \frac{\phi}{2}\right)")
            st.latex(r"\sigma_{total} = \sigma_{soil} + \sigma_q + \sigma_w")
            st.markdown(f"- **Friction Angle ($\phi$):** {phi}° &nbsp; | &nbsp; **Ka:** {Ka:.4f}")
            st.markdown(f"- **Total lateral resultant:** {total_lateral_resultant:.3f} kN")
            st.markdown(f"- **Total soil vertical (on heel):** {soil_total:.3f} kN")
            
            st.markdown("---")
            st.markdown("#### Concrete Self-Weight")
            st.latex(r"W_{conc} = \gamma_c\,A_{conc}\,t")
            st.markdown(f"- **Concrete Weight (W_conc):** {W_conc:.3f} kN")
            st.markdown(f"- **Stem Area:** {wall_area:.3f} m², &nbsp; **Footing Area:** {footing_area:.3f} m²")
            
            st.markdown("---")
            st.markdown("#### Water Uplift Pressure")
            st.latex(r"F_{uplift,node} = (\gamma_w H_{wtr})\cdot L\cdot t")
            st.markdown(f"- **Water uplift force:** {uplift_total:.3f} kN")
            st.markdown(f"- **Net vertical load:** {net_vertical:.3f} kN")
            
            st.markdown("---")
            st.markdown("#### Seismic Loading (Pseudo-static)")
            st.latex(r"F_{seismic} = k_h\,(W_{soil} + W_{conc})")
            st.markdown(f"- **Seismic Coefficient ($k_h$):** {kh:.4f}")
            st.markdown(f"- **Pseudo-static lateral force:** {F_seismic_pseudo:.3f} kN")
            st.markdown(f"- **Transient Peak Acceleration estimate:** {peak_accel:.3f} m/s²")
            st.markdown(f"- **Inertial Peak Force estimate:** {inertial_peak:.3f} kN")


# --- RIGHT COLUMN: ANALYSIS RESULTS ---
with col_right:

    # Card C: Stem Bending Moment Diagram
    with st.container(border=True):
        st.subheader("📉 Stem Bending Moment Diagram")
        st.markdown('<div class="section-desc">Internal bending moment contour along the height of the stem wall.</div>', unsafe_allow_html=True)
        
        try:
            node_moments_dict = get_stem_moments()
            nds_val_moment = np.zeros(total_nodes)
            for i in range(1, total_nodes + 1):
                if i in node_moments_dict:
                    nds_val_moment[i-1] = node_moments_dict[i]
                else:
                    nds_val_moment[i-1] = 0.0
            
            plt.close('all')
            opsv.plot_stress_2d(nds_val_moment)
            fig_m = plt.gcf()
            fig_m.set_size_inches(6, 4.5)
            try:
                ax = fig_m.axes[0]
                ax.set_xlabel('x [m]')
                ax.set_ylabel('y [m]')
                ax.set_title('Stem Node Moment (kN·m)')
                ax.relim()
                ax.autoscale_view()
            except Exception:
                pass
            fig_m.tight_layout()
            st.pyplot(fig_m)
            plt.close(fig_m)
        except Exception as e:
            st.error(f'Failed to compute/plot stem node moments: {e}')

    # Card D: Finite Element Stress Contours
    with st.container(border=True):
        st.subheader("🎭 Finite Element Stress Contours")
        st.markdown('<div class="section-desc">2D stress contour plots mapped across the retaining wall mesh elements.</div>', unsafe_allow_html=True)
        
        stress_col1, stress_col2 = st.columns(2)
        
        with stress_col1:
            st.markdown("**X-Direction stress ($\sigma_{xx}$)**")
            plt.close('all')
            try:
                sig_out = opsv.sig_out_per_node()
                j, jstr = 0, 'σ_xx (kN/m^2)'
                nds_val = np.array(sig_out)[:, j]

                opsv.plot_stress_2d(nds_val)
                fig_res = plt.gcf()
                fig_res.set_size_inches(4.5, 3.8)
                try:
                    ax = fig_res.axes[0]
                    ax.set_xlabel('x [m]')
                    ax.set_ylabel('y [m]')
                    ax.relim()
                    ax.autoscale_view()
                except Exception:
                    pass
                fig_res.tight_layout()
                st.pyplot(fig_res)
                plt.close(fig_res)
            except Exception as e:
                st.error(f'Failed to render σ_xx: {e}')
                
        with stress_col2:
            st.markdown("**Y-Direction stress ($\sigma_{yy}$)**")
            plt.close('all')
            try:
                sig_out = opsv.sig_out_per_node()
                j, jstr = 1, 'σ_yy (kN/m^2)'
                nds_val = np.array(sig_out)[:, j]

                opsv.plot_stress_2d(nds_val)
                fig_res2 = plt.gcf()
                fig_res2.set_size_inches(4.5, 3.8)
                try:
                    ax2 = fig_res2.axes[0]
                    ax2.set_xlabel('x [m]')
                    ax2.set_ylabel('y [m]')
                    ax2.relim()
                    ax2.autoscale_view()
                except Exception:
                    pass
                fig_res2.tight_layout()
                st.pyplot(fig_res2)
                plt.close(fig_res2)
            except Exception as e:
                st.error(f'Failed to render σ_yy: {e}')

    # Card E: Detailed Stability Verification
    with st.container(border=True):
        st.subheader("🛡️ Stability Verification Report")
        st.markdown('<div class="section-desc">Detailed overturning, sliding, and stem junction stability verification checklist.</div>', unsafe_allow_html=True)
        
        ot_limit = 1.5
        slide_limit = 1.5
        
        def status_badge_html(val, limit):
            if val is None:
                return "<span style='color:gray; font-weight:bold;'>N/A</span>"
            if val == float('inf'):
                return "<span style='background-color:#28a745; color:white; padding:3px 10px; border-radius:12px; font-weight:bold; font-size:12px; text-transform:uppercase;'>PASS (∞)</span>"
            if val >= limit:
                return f"<span style='background-color:#28a745; color:white; padding:3px 10px; border-radius:12px; font-weight:bold; font-size:12px; text-transform:uppercase;'>PASS ({val:.3f})</span>"
            else:
                return f"<span style='background-color:#dc3545; color:white; padding:3px 10px; border-radius:12px; font-weight:bold; font-size:12px; text-transform:uppercase;'>FAIL ({val:.3f})</span>"

        fs_ot_global_val = f"{FS_ot_global:.3f}" if FS_ot_global is not None else "0.000"
        fs_slide_global_val = f"{FS_slide_global:.3f}" if FS_slide_global is not None else "0.000"
        fs_ot_stem_val = f"{FS_ot_stem:.3f}" if FS_ot_stem is not None else "0.000"

        report_md = f"""
        - **Global Overturning Safety** (pivot $x=0.0$):
          - $M_{{driving}}$ = **{M_ot_global:.3f} kN·m** &nbsp; | &nbsp; $M_{{resisting}}$ = **{M_res_global:.3f} kN·m**
          - **FS = {fs_ot_global_val}** &nbsp;&nbsp; {status_badge_html(FS_ot_global, ot_limit)}
        
        - **Global Sliding Safety**:
          - $R_{{sliding}}$ = **{R_slide:.3f} kN** &nbsp; | &nbsp; $F_{{driving}}$ = **{total_lateral_resultant:.3f} kN**
          - **FS = {fs_slide_global_val}** &nbsp;&nbsp; {status_badge_html(FS_slide_global, slide_limit)}
          
        - **Stem Junction Overturning Safety** (pivot $x={x_pivot_stem:.3f}$):
          - $M_{{driving}}$ = **{M_ot_stem:.3f} kN·m** &nbsp; | &nbsp; $M_{{resisting}}$ = **{M_res_stem:.3f} kN·m**
          - **FS = {fs_ot_stem_val}** &nbsp;&nbsp; {status_badge_html(FS_ot_stem, ot_limit)}
        """
        if enable_bearing:
            fs_bearing_val = f"{FS_actual_bearing:.3f}" if FS_actual_bearing != float('inf') else "∞"
            report_md += f"""
        - **Bearing Capacity Safety**:
          - $\\sigma_{{max}}$ = **{sigma_max:.3f} kPa** &nbsp; | &nbsp; $q_{{allow}}$ = **{q_allow:.3f} kPa** (SF = {FS_bearing:.2f})
          - **FS = {fs_bearing_val}** &nbsp;&nbsp; {status_badge_html(FS_actual_bearing, FS_bearing)}
        """
        st.markdown(report_md, unsafe_allow_html=True)
        
        with st.expander("📖 View Overturning and Sliding Equations"):
            st.markdown("**Global Overturning Verification**")
            st.latex(r"M_{ot} = \sum |F_{x,i} \cdot z_i|")
            st.latex(r"M_{res} = W_{base}\cdot x_{base} + W_{stem,rect}\cdot x_{stem,rect} + W_{stem,tri}\cdot x_{stem,tri} + W_{soil}\cdot x_{soil} - M_{uplift}")
            st.latex(r"FS_{ot} = \frac{M_{res}}{M_{ot}} \ge 1.5")
            
            st.markdown("**Global Sliding Verification**")
            st.latex(r"R_{slide} = W_{down} \cdot \tan(\phi)")
            st.latex(r"FS_{slide} = \frac{R_{slide}}{F_{drive}} \ge 1.5")
            
            st.markdown("**Stem Overturning Verification**")
            st.latex(r"M_{ot} = \sum |F_{x,i} \cdot (x_i - x_{pivot})|")
            st.latex(r"M_{res} = W_{down} \cdot |x_{centroid} - x_{pivot}|")
            st.latex(r"FS_{ot} = \frac{M_{res}}{M_{ot}} \ge 1.5")

            if enable_bearing:
                st.markdown("**Bearing Capacity Verification**")
                st.latex(r"\sigma_{max} = \max(q_{toe}, q_{heel})")
                st.latex(r"q_{allow} = \frac{q_{ult}}{FS_{bearing}}")
                st.latex(r"FS_{actual} = \frac{q_{ult}}{\sigma_{max}} \ge FS_{bearing}")

    # Card F: Concrete Strength Verification
    with st.container(border=True):
        st.subheader("🛡️ Concrete Strength Verification")
        st.markdown('<div class="section-desc">Comparison of structural bending moment demands against ultimate cross-section capacities.</div>', unsafe_allow_html=True)
        
        def check_badge(is_pass):
            if is_pass:
                return "<span style='background-color:#28a745; color:white; padding:3px 10px; border-radius:12px; font-weight:bold; font-size:12px;'>PASS</span>"
            else:
                return "<span style='background-color:#dc3545; color:white; padding:3px 10px; border-radius:12px; font-weight:bold; font-size:12px;'>FAIL</span>"
                
        # Comparison Table using standard HTML
        st.markdown(f"""
        <table style="width:100%; border-collapse: collapse; margin-bottom: 15px;">
            <thead>
                <tr style="border-bottom: 2px solid #eef1f6; text-align: left;">
                    <th style="padding: 8px;">Component</th>
                    <th style="padding: 8px;">Bending Type</th>
                    <th style="padding: 8px; text-align: right;">Actual M* (kN·m)</th>
                    <th style="padding: 8px; text-align: right;">Capacity φMu (kN·m)</th>
                    <th style="padding: 8px; text-align: center;">Status</th>
                </tr>
            </thead>
            <tbody>
                <tr style="border-bottom: 1px solid #eef1f6;">
                    <td style="padding: 8px; font-weight: 600;">Stem Wall (Base)</td>
                    <td style="padding: 8px; color: #6c757d; font-size: 0.9rem;">Bending major</td>
                    <td style="padding: 8px; text-align: right; font-family: monospace;">{M_actual_stem:.2f}</td>
                    <td style="padding: 8px; text-align: right; font-family: monospace;">{M_u_stem_allow:.2f}</td>
                    <td style="padding: 8px; text-align: center;">{check_badge(stem_pass)}</td>
                </tr>
                <tr style="border-bottom: 1px solid #eef1f6;">
                    <td style="padding: 8px; font-weight: 600;">Footing Toe (Slab)</td>
                    <td style="padding: 8px; color: #6c757d; font-size: 0.9rem;">Sagging (Tension Bottom)</td>
                    <td style="padding: 8px; text-align: right; font-family: monospace;">{M_actual_toe:.2f}</td>
                    <td style="padding: 8px; text-align: right; font-family: monospace;">{M_u_ftg_sag_allow:.2f}</td>
                    <td style="padding: 8px; text-align: center;">{check_badge(toe_pass)}</td>
                </tr>
                <tr style="border-bottom: 1px solid #eef1f6;">
                    <td style="padding: 8px; font-weight: 600;">Footing Heel (Slab)</td>
                    <td style="padding: 8px; color: #6c757d; font-size: 0.9rem;">Hogging (Tension Top)</td>
                    <td style="padding: 8px; text-align: right; font-family: monospace;">{M_actual_heel:.2f}</td>
                    <td style="padding: 8px; text-align: right; font-family: monospace;">{M_u_ftg_hog_allow:.2f}</td>
                    <td style="padding: 8px; text-align: center;">{check_badge(heel_pass)}</td>
                </tr>
            </tbody>
        </table>
        """, unsafe_allow_html=True)
        
        # Conclusion box
        all_concrete_pass = stem_pass and toe_pass and heel_pass
        conclusion_info = f" (using Main Rebar **D{stem_rebar_dia}** for Stem Wall, and **D{ftg_rebar_dia}** for Footing Slabs)"
        if all_concrete_pass:
            st.success(f"🤖 **Conclusion**: The structural verification is successful! All key components (Stem Wall, Footing Toe, and Footing Heel) satisfy the bending capacity checks. The current reinforcement configuration{conclusion_info} is **adequate** and **safe** under the design limit state load combinations.")
        else:
            failed_components = []
            if not stem_pass: failed_components.append("Stem Wall")
            if not toe_pass: failed_components.append("Footing Toe")
            if not heel_pass: failed_components.append("Footing Heel")
            st.error(f"⚠️ **Conclusion**: Structural verification **FAILED** due to insufficient bending capacity in: **{', '.join(failed_components)}**. The current reinforcement configuration{conclusion_info} is **inadequate**. Please consider increasing the rebar diameter, decreasing rebar spacing, or increasing the thickness of the failing member.")

        # Expander to see layouts
        with st.expander("🔍 View Concrete Rebar Cross Sections"):
            if fig_stem_sec is not None:
                st.pyplot(fig_stem_sec)
            if fig_ftg_sec is not None:
                st.pyplot(fig_ftg_sec)

# ----------------------------------------------------
# Pile Foundation Analysis Results Section
# ----------------------------------------------------
if enable_pile:
    st.markdown("---")
    st.markdown("<h2 style='text-align: center; color: #1f4068;'>🧱 Pile Foundation Analysis & Design</h2>", unsafe_allow_html=True)
    st.markdown('<div class="section-desc" style="text-align: center; margin-bottom: 24px;">Comprehensive evaluation of base anchor piles under combined vertical load, lateral shear, and bending demand.</div>', unsafe_allow_html=True)
    
    pile_kpi_cols = st.columns(6)
    with pile_kpi_cols[0]:
        with st.container(border=True):
            st.metric("Toe Pile Axial Force", f"{P_toe:.2f} kN", delta="Compression" if P_toe >= 0 else "Tension", delta_color="normal" if P_toe >= 0 else "inverse")
    with pile_kpi_cols[1]:
        with st.container(border=True):
            st.metric("Toe Pile Max Shear", f"{V_max_toe:.2f} kN")
    with pile_kpi_cols[2]:
        with st.container(border=True):
            st.metric("Toe Pile Max Moment", f"{M_max_toe:.2f} kN·m")
    with pile_kpi_cols[3]:
        with st.container(border=True):
            st.metric("Heel Pile Axial Force", f"{P_heel:.2f} kN", delta="Compression" if P_heel >= 0 else "Tension", delta_color="normal" if P_heel >= 0 else "inverse")
    with pile_kpi_cols[4]:
        with st.container(border=True):
            st.metric("Heel Pile Max Shear", f"{V_max_heel:.2f} kN")
    with pile_kpi_cols[5]:
        with st.container(border=True):
            st.metric("Heel Pile Max Moment", f"{M_max_heel:.2f} kN·m")

    # Layout for plots
    plot_col1, plot_col2, plot_col3 = st.columns([1.2, 1, 1])
    
    with plot_col1:
        with st.container(border=True):
            st.subheader("📐 Pile Foundation Drawing")
            st.markdown('<div class="section-desc">Retaining wall section showing embedded base anchors with reaction forces.</div>', unsafe_allow_html=True)
            # Create SVG drawing of piles below footing base
            try:
                mult = 1000
                max_h = max(Hw, h_soil, h_soil_toe, Hwtr, Hwtr_front) + h_ftg
                # Set display parameters to shrink canvas borders and margins
                L_draw = 1.2
                margin_bottom = 0.4
                margin_top = 0.4
                margin_left = 0.4
                margin_right = 0.4
                
                W_canvas = ftg + margin_left + margin_right
                H_canvas = max_h + L_draw + margin_bottom + margin_top
                
                d_pile_svg = draw.Drawing(W_canvas * mult, H_canvas * mult, origin='bottom-left')
                
                x_shift = margin_left * mult
                yshift_pile = (margin_bottom + L_draw) * mult
                
                # Draw soil - Dry soil (above water table)
                if h_soil > Hwtr:
                    d_pile_svg.append(draw.Lines(x_shift + (toe+bot_wall)*mult, -(h_wet_height+h_ftg)*mult - yshift_pile,
                                            x_shift + ftg*mult, -(h_wet_height+h_ftg)*mult - yshift_pile,
                                            x_shift + ftg*mult, -(h_soil+h_ftg)*mult - yshift_pile,
                                            x_shift + (toe+bot_wall)*mult, -(h_soil+h_ftg)*mult - yshift_pile,
                                            close=True, fill='#5BC2A5', stroke='black', stroke_width=3))

                if h_wet_height > 0:
                    d_pile_svg.append(draw.Lines(x_shift + ftg*mult, -h_ftg*mult - yshift_pile,
                                            x_shift + ftg*mult, -(h_wet_height+h_ftg)*mult - yshift_pile,
                                            x_shift + (toe+bot_wall)*mult, -(h_wet_height+h_ftg)*mult - yshift_pile,
                                            x_shift + (toe+bot_wall)*mult, -h_ftg*mult - yshift_pile,
                                            close=True, fill='#A6F527', stroke='black', stroke_width=3))

                if h_soil_toe > 0:
                    d_pile_svg.append(draw.Lines(x_shift + 0, -h_ftg*mult - yshift_pile,
                                            x_shift + toe*mult, -h_ftg*mult - yshift_pile,
                                            x_shift + toe*mult, -(h_ftg + h_soil_toe)*mult - yshift_pile,
                                            x_shift + 0, -(h_ftg + h_soil_toe)*mult - yshift_pile,
                                            close=True, fill='#5BC2A5', stroke='black', stroke_width=3))

                # Draw Retaining wall (shifted up by L_draw + margin_bottom)
                d_pile_svg.append(draw.Lines(x_shift + 0.0,  -yshift_pile,
                                        x_shift + 0.0,  -h_ftg*mult - yshift_pile,
                                        x_shift + toe*mult, -h_ftg*mult - yshift_pile,
                                        x_shift + (toe+(taper/4))*mult, -(h_ftg+(Hw/4))*mult - yshift_pile,
                                        x_shift + (toe+(taper/2))*mult, -(h_ftg+(Hw/2))*mult - yshift_pile,
                                        x_shift + (toe+(taper*(3/4)))*mult, -(h_ftg+(Hw*(3/4)))*mult - yshift_pile,
                                        x_shift + (toe+taper)*mult, -(h_ftg+Hw)*mult - yshift_pile,
                                        x_shift + (toe+bot_wall)*mult, -(h_ftg+Hw)*mult - yshift_pile,
                                        x_shift + (toe+bot_wall)*mult, -(h_ftg+(Hw*(3/4)))*mult - yshift_pile,
                                        x_shift + (toe+bot_wall)*mult, -(h_ftg+(Hw/2))*mult - yshift_pile,
                                        x_shift + (toe+bot_wall)*mult, -(h_ftg+(Hw/4))*mult - yshift_pile,
                                        x_shift + (toe+bot_wall)*mult, -h_ftg*mult - yshift_pile,
                                        x_shift + ftg*mult, -h_ftg*mult - yshift_pile,
                                        x_shift + ftg*mult, -yshift_pile,
                                        close=True, fill='#eeee00', stroke='black', stroke_width=50))
                
                # Draw Piles with zig-zag break lines
                w_draw = (diameter_pile if pile_shape == 'Circle' else width_x_pile)
                y_top = -yshift_pile
                y_bot = -margin_bottom * mult
                zag1_y = y_bot + 0.05 * mult
                zag2_y = y_bot - 0.05 * mult
                
                # Toe Pile
                x_center_toe = x_shift + pile_offset * mult
                x_l_toe = x_center_toe - (w_draw / 2.0) * mult
                x_r_toe = x_center_toe + (w_draw / 2.0) * mult
                d_pile_svg.append(draw.Lines(
                    x_l_toe, y_top,
                    x_r_toe, y_top,
                    x_r_toe, y_bot,
                    x_l_toe + (w_draw * 0.75) * mult, zag1_y,
                    x_l_toe + (w_draw * 0.25) * mult, zag2_y,
                    x_l_toe, y_bot,
                    close=True,
                    fill='#CCCCCC' if pile_material == 'Concrete' else '#777777', 
                    stroke='black', stroke_width=20
                ))
                
                # Heel Pile
                x_center_heel = x_shift + (ftg - pile_offset) * mult
                x_l_heel = x_center_heel - (w_draw / 2.0) * mult
                x_r_heel = x_center_heel + (w_draw / 2.0) * mult
                d_pile_svg.append(draw.Lines(
                    x_l_heel, y_top,
                    x_r_heel, y_top,
                    x_r_heel, y_bot,
                    x_l_heel + (w_draw * 0.75) * mult, zag1_y,
                    x_l_heel + (w_draw * 0.25) * mult, zag2_y,
                    x_l_heel, y_bot,
                    close=True,
                    fill='#CCCCCC' if pile_material == 'Concrete' else '#777777', 
                    stroke='black', stroke_width=20
                ))
                
                # Draw pile reaction force arrows at the top of the piles
                def draw_reaction_arrow(d_draw, x_c, p_val):
                    color = 'red' if p_val >= 0 else 'blue'
                    y_top = -yshift_pile
                    if p_val >= 0:
                        d_draw.append(draw.Line(x_c, y_top + 300, x_c, y_top + 50, stroke=color, stroke_width=15))
                        d_draw.append(draw.Lines(x_c - 50, y_top + 150, x_c, y_top + 50, x_c + 50, y_top + 150, stroke=color, stroke_width=15))
                    else:
                        d_draw.append(draw.Line(x_c, y_top + 50, x_c, y_top + 300, stroke=color, stroke_width=15))
                        d_draw.append(draw.Lines(x_c - 50, y_top + 200, x_c, y_top + 300, x_c + 50, y_top + 200, stroke=color, stroke_width=15))
                    d_draw.append(draw.Text(f"{abs(p_val):.1f} kN", 180, x_c + 120, y_top + 220, fill='black', font_weight='bold'))
                
                draw_reaction_arrow(d_pile_svg, x_center_toe, P_toe)
                draw_reaction_arrow(d_pile_svg, x_center_heel, P_heel)
                
                # Labels
                d_pile_svg.append(draw.Text("TOE", 180, x_center_toe, y_top + 600, fill='black', text_anchor='middle', font_weight='bold'))
                d_pile_svg.append(draw.Text("HEEL", 180, x_center_heel, y_top + 600, fill='black', text_anchor='middle', font_weight='bold'))
                
                d_pile_svg.set_render_size(600, 500)
                st.image(d_pile_svg.as_svg(), use_container_width=True)
            except Exception as e_svg:
                st.error(f"Error drawing piles SVG: {e_svg}")
                
    with plot_col2:
        with st.container(border=True):
            st.subheader("📈 Winkler Lateral Response")
            st.markdown('<div class="section-desc">Pile deflection, shear force, and bending moment profiles.</div>', unsafe_allow_html=True)
            try:
                fig_winkler, (ax_d, ax_v, ax_m) = plt.subplots(1, 3, figsize=(10, 8), sharey=True)
                
                elevs_defl_toe = res_toe.deflection['Elevation [m]']
                defl_toe = res_toe.deflection['Deflection [m]'] * 1000.0
                elevs_forces_toe = res_toe.forces['Elevation [m]']
                shear_toe = res_toe.forces['V [kN]']
                moment_toe = res_toe.forces['M [kNm]']
                
                elevs_defl_heel = res_heel.deflection['Elevation [m]']
                defl_heel = res_heel.deflection['Deflection [m]'] * 1000.0
                elevs_forces_heel = res_heel.forces['Elevation [m]']
                shear_heel = res_heel.forces['V [kN]']
                moment_heel = res_heel.forces['M [kNm]']
                
                ax_d.plot(defl_toe, elevs_defl_toe, 'b-', label='Toe Pile')
                ax_d.plot(defl_heel, elevs_defl_heel, 'r--', label='Heel Pile')
                ax_d.set_xlabel('Deflection (mm)')
                ax_d.set_ylabel('Elevation (m)')
                ax_d.grid(True)
                ax_d.legend(loc='lower right')
                
                ax_v.plot(shear_toe, elevs_forces_toe, 'b-', label='Toe Pile')
                ax_v.plot(shear_heel, elevs_forces_heel, 'r--', label='Heel Pile')
                ax_v.set_xlabel('Shear Force (kN)')
                ax_v.grid(True)
                
                ax_m.plot(moment_toe, elevs_forces_toe, 'b-', label='Toe Pile')
                ax_m.plot(moment_heel, elevs_forces_heel, 'r--', label='Heel Pile')
                ax_m.set_xlabel('Moment (kN·m)')
                ax_m.grid(True)
                
                fig_winkler.tight_layout()
                st.pyplot(fig_winkler)
            except Exception as e_w:
                st.error(f"Error plotting Winkler response: {e_w}")
                
    with plot_col3:
        with st.container(border=True):
            st.subheader("🎯 Moment-Axial Interaction")
            st.markdown('<div class="section-desc">Pile design capacity envelope showing the demand point for each pile.</div>', unsafe_allow_html=True)
            try:
                fig_int, ax_int = plt.subplots(figsize=(6, 5))
                ax_int.plot(curve_m, curve_n, 'k-', linewidth=2, label='Capacity Envelope')
                ax_int.scatter([M_max_toe], [P_toe], color='blue', marker='o', s=100, label=f'Toe Pile (M={M_max_toe:.1f}, P={P_toe:.1f})')
                ax_int.scatter([M_max_heel], [P_heel], color='red', marker='^', s=100, label=f'Heel Pile (M={M_max_heel:.1f}, P={P_heel:.1f})')
                ax_int.axhline(0, color='gray', linestyle='--', linewidth=0.8)
                ax_int.axvline(0, color='gray', linestyle='--', linewidth=0.8)
                ax_int.set_xlabel('Bending Moment (kN·m)')
                ax_int.set_ylabel('Axial Force (kN)')
                ax_int.set_title('Pile Bending-Axial Interaction')
                ax_int.legend(loc='upper right')
                ax_int.grid(True)
                fig_int.tight_layout()
                st.pyplot(fig_int)
                plt.close(fig_int)
            except Exception as e_i:
                st.error(f"Error plotting Interaction Diagram: {e_i}")

    # Pile Cross Section Geometry
    st.markdown("---")
    st.subheader("🔍 Pile Cross Section Geometry")
    if fig_pile_sec is not None:
        col_p1, col_p2 = st.columns([1, 1])
        with col_p1:
            st.pyplot(fig_pile_sec)
        with col_p2:
            if pile_material == 'Concrete':
                if pile_shape == 'Circle':
                    st.markdown(f"""
                    **Circular Concrete Pile Section**:
                    - **Pile Diameter**: `{diameter_pile:.3f} m` ({diameter_pile*1000:.0f} mm)
                    - **Concrete Cover**: `{cover_pile:.1f} mm`
                    - **Reinforcement**: **{int(n_rebar_pile)}** main bars of **D{rebar_dia_pile}**
                    - **Friction Array Radius**: `{diameter_pile*1000/2 - cover_pile:.1f} mm`
                    """)
                else:
                    st.markdown(f"""
                    **Rectangular Concrete Pile Section**:
                    - **Pile Dimensions (dx × dy)**: `{width_x_pile:.3f} m × {width_y_pile:.3f} m` ({width_x_pile*1000:.0f} mm × {width_y_pile*1000:.0f} mm)
                    - **Concrete Cover**: `{cover_pile:.1f} mm`
                    - **Reinforcement**: **4** corner bars of **D{rebar_dia_pile}**
                    """)
            else:
                if pile_shape == 'Circle':
                    st.markdown(f"""
                    **Solid Steel Circular Pile Section**:
                    - **Pile Diameter**: `{diameter_pile:.3f} m` ({diameter_pile*1000:.0f} mm)
                    """)
                else:
                    st.markdown(f"""
                    **Solid Steel Rectangular Pile Section**:
                    - **Pile Dimensions (dx × dy)**: `{width_x_pile:.3f} m × {width_y_pile:.3f} m` ({width_x_pile*1000:.0f} mm × {width_y_pile*1000:.0f} mm)
                    """)

    # Pile Foundation Report Cards
    st.subheader("📋 Pile Foundation Verification Report")
    
    def status_span(ok):
        if ok:
            return "<span style='background-color:#28a745; color:white; padding:3px 10px; border-radius:12px; font-weight:bold; font-size:12px;'>PASS</span>"
        return "<span style='background-color:#dc3545; color:white; padding:3px 10px; border-radius:12px; font-weight:bold; font-size:12px;'>FAIL</span>"
        
    st.markdown(f"""
    <table style="width:100%; border-collapse: collapse; margin-bottom: 25px;">
        <thead>
            <tr style="border-bottom: 2px solid #eef1f6; text-align: left;">
                <th style="padding: 10px;">Check Item</th>
                <th style="padding: 10px;">Demand (Actual)</th>
                <th style="padding: 10px;">Capacity (Allowable)</th>
                <th style="padding: 10px; text-align: center;">Safety Factor</th>
                <th style="padding: 10px; text-align: center;">Status</th>
            </tr>
        </thead>
        <tbody>
            <tr style="border-bottom: 1px solid #eef1f6;">
                <td style="padding: 10px; font-weight: bold;">Toe Pile Axial (Compression)</td>
                <td style="padding: 10px; font-family: monospace;">{P_toe:.2f} kN</td>
                <td style="padding: 10px; font-family: monospace;">{Q_comp_allow:.2f} kN</td>
                <td style="padding: 10px; text-align: center; font-family: monospace;">FS = {FS_pile_axial:.2f}</td>
                <td style="padding: 10px; text-align: center;">{status_span(pile_pass_toe_axial)}</td>
            </tr>
            <tr style="border-bottom: 1px solid #eef1f6;">
                <td style="padding: 10px; font-weight: bold;">Toe Pile Lateral Shear</td>
                <td style="padding: 10px; font-family: monospace;">{V_max_toe:.2f} kN</td>
                <td style="padding: 10px; font-family: monospace;">{H_allow:.2f} kN</td>
                <td style="padding: 10px; text-align: center; font-family: monospace;">FS = {FS_pile_lateral:.2f}</td>
                <td style="padding: 10px; text-align: center;">{status_span(pile_pass_toe_lateral)}</td>
            </tr>
            <tr style="border-bottom: 1px solid #eef1f6;">
                <td style="padding: 10px; font-weight: bold;">Toe Pile Structural (Bending-Axial)</td>
                <td style="padding: 10px; font-family: monospace;">M={M_max_toe:.1f} kNm, P={P_toe:.1f} kN</td>
                <td style="padding: 10px; font-family: monospace;">Interaction Envelope</td>
                <td style="padding: 10px; text-align: center; font-family: monospace;">-</td>
                <td style="padding: 10px; text-align: center;">{status_span(pile_pass_toe_interaction)}</td>
            </tr>
            <tr style="border-bottom: 1px solid #eef1f6;">
                <td style="padding: 10px; font-weight: bold;">Heel Pile Axial ({"Compression" if P_heel >= 0 else "Tension"})</td>
                <td style="padding: 10px; font-family: monospace;">{abs(P_heel):.2f} kN</td>
                <td style="padding: 10px; font-family: monospace;">{Q_comp_allow if P_heel >= 0 else Q_tens_allow:.2f} kN</td>
                <td style="padding: 10px; text-align: center; font-family: monospace;">FS = {FS_pile_axial:.2f}</td>
                <td style="padding: 10px; text-align: center;">{status_span(pile_pass_heel_axial)}</td>
            </tr>
            <tr style="border-bottom: 1px solid #eef1f6;">
                <td style="padding: 10px; font-weight: bold;">Heel Pile Lateral Shear</td>
                <td style="padding: 10px; font-family: monospace;">{V_max_heel:.2f} kN</td>
                <td style="padding: 10px; font-family: monospace;">{H_allow:.2f} kN</td>
                <td style="padding: 10px; text-align: center; font-family: monospace;">FS = {FS_pile_lateral:.2f}</td>
                <td style="padding: 10px; text-align: center;">{status_span(pile_pass_heel_lateral)}</td>
            </tr>
            <tr style="border-bottom: 1px solid #eef1f6;">
                <td style="padding: 10px; font-weight: bold;">Heel Pile Structural (Bending-Axial)</td>
                <td style="padding: 10px; font-family: monospace;">M={M_max_heel:.1f} kNm, P={P_heel:.1f} kN</td>
                <td style="padding: 10px; font-family: monospace;">Interaction Envelope</td>
                <td style="padding: 10px; text-align: center; font-family: monospace;">-</td>
                <td style="padding: 10px; text-align: center;">{status_span(pile_pass_heel_interaction)}</td>
            </tr>
        </tbody>
    </table>
    """, unsafe_allow_html=True)
