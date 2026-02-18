# CALCULATION OF THE ANGLES AND THE PAIRS
import numpy as np

# =====================
# DETECTOR GEOMETRY
# =====================
a = 8.2  # Side of the squared front face of the clover in cm
b = 7.0  # Depth in cm of the HpGe crystals
"""
# Angular deviation from the clover axis to origin of coordinates
tilt = np.array([0.377, 0.055, 0.386, 0.147, 0.255, 0.728, 0.396, 0.033, 0.746])

# Rotation in the ZY plane of the clover
roll_angle = np.array([0.0193, 0.0178, -0.2538, 0.0152, -0.0577, 0.0200, -0.0125, -0.0013, -0.3315])

# End vectors indexed by clover
end_vectors = [
    np.array([ -0.170, -46.590, -46.130]),
    np.array([ -0.310, -84.320,  -0.040]),
    np.array([ 0.120, -60.760,  60.370]),
    np.array([ 63.400,   0.220,  63.180]),
    np.array([ 59.440,   0.200,   0.020]),
    np.array([ 60.410,  -0.110, -59.440]),
    np.array([ -0.380,  45.460,  44.960]),
    np.array([ -0.010,  84.120,  -0.060]),
    np.array([ -0.110,  43.440, -42.830])
]

# Theta, phi and distance from implantation point to center of the face of each clover
clovers = {
    'g2_A': (np.radians(135.11), np.radians(269.7962), 12.7989),
    'g2_B': (np.radians(90.0037), np.radians(269.839), 12.4805),
    'g2_C': (np.radians(44.9839), np.radians(270.5336), 12.7650),
    'g3_A': (np.radians(45.0469), np.radians(0.0024), 12.6100),
    'g3_B': (np.radians(89.9437), np.radians(359.9342), 12.5027),
    'g3_C': (np.radians(135.0395), np.radians(359.1597), 13.3711),
    'g4_A': (np.radians(44.9311), np.radians(90.5476), 12.8691),
    'g4_B': (np.radians(90.072), np.radians(89.6847), 12.5464),
    'g4_C': (np.radians(135.0007), np.radians(89.2739), 12.8949)
}
"""

tilt = np.array([1.695, 1.075, 0.377, 0.055, 0.386, 0.147, 0.255, 0.728, 0.396, 0.033, 0.746, 1.529, 1.587])

roll_angle = np.array([-0.2885, 1.0288, 0.0193, 0.0178, -0.2538, 0.0152, -0.0577, 0.0200, -0.0125, -0.0013, -0.3315, 0.4523, -0.5372])

# Direcciones de referencia para aplicar el tilt
end_vectors = [
    np.array([-49.470, -36.400,  61.410]),
    np.array([-50.890, -36.900, -60.400]),
    np.array([ -0.170, -46.590, -46.130]),
    np.array([ -0.310, -84.320,  -0.040]),
    np.array([ 0.120, -60.760,  60.370]),
    np.array([ 63.400,   0.220,  63.180]),
    np.array([ 59.440,   0.200,   0.020]),
    np.array([ 60.410,  -0.110, -59.440]),
    np.array([ -0.380,  45.460,  44.960]),
    np.array([ -0.010,  84.120,  -0.060]),
    np.array([ -0.110,  43.440, -42.830]),
    np.array([-38.900,  27.840,  46.040]),
    np.array([-39.530,  28.700, -46.540])
]

# =====================
# POSICIONES DE LOS CLOVERS (θ, φ, r)
# =====================
clovers = {
    'g1_A': (np.radians(43.3236), np.radians(216.8214), 15.6368),
    'g1_B': (np.radians(134.8903), np.radians(215.5797), 15.7215),
    'g2_A': (np.radians(135.11), np.radians(269.7962), 12.7989),
    'g2_B': (np.radians(90.0037), np.radians(269.839), 12.4805),
    'g2_C': (np.radians(44.9839), np.radians(270.5336), 12.7650),
    'g3_A': (np.radians(45.0469), np.radians(0.0024), 12.6100),
    'g3_B': (np.radians(89.9437), np.radians(359.9342), 12.5027),
    'g3_C': (np.radians(135.0395), np.radians(359.1597), 13.3711),
    'g4_A': (np.radians(44.9311), np.radians(90.5476), 12.8691),
    'g4_B': (np.radians(90.072), np.radians(89.6847), 12.5464),
    'g4_C': (np.radians(135.0007), np.radians(89.2739), 12.8949),
    'g5_A': (np.radians(44.5706), np.radians(144.6118), 15.0870),
    'g5_B': (np.radians(135.1691), np.radians(143.5214), 15.5307)
}

# Offsets de cristales en sistema local del clover
crystal_offsets = {
    'yellow': np.array([0, -a / 4, -a / 4]),
    'green':  np.array([0,  a / 4, -a / 4]),
    'blue':   np.array([0,  a / 4,  a / 4]),
    'red':    np.array([0, -a / 4,  a / 4]),
}
original_colors = ['yellow', 'green', 'blue', 'red']

# =====================
# FUNCIONES AUXILIARES
# =====================
def spherical_to_cartesian(theta, phi, r, depth=b):
    R = r + depth / 2
    x = R * np.sin(theta) * np.cos(phi)
    y = R * np.sin(theta) * np.sin(phi)
    z = R * np.cos(theta)
    return np.array([x, y, z])

def rotation_matrix(axis, alpha):
    axis = axis / np.linalg.norm(axis)
    a = np.cos(alpha/2)
    b, c, d = -axis * np.sin(alpha/2)
    return np.array([
        [a*a + b*b - c*c - d*d, 2*(b*c - a*d),     2*(b*d + a*c)],
        [2*(b*c + a*d),         a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
        [2*(b*d - a*c),         2*(c*d + a*b),     a*a + d*d - b*b - c*c]
    ])

def build_local_basis_face(face_normal):
    x_local = face_normal / np.linalg.norm(face_normal)
    y_local = np.cross([0,0,1], x_local)
    if np.linalg.norm(y_local) < 1e-6:
        y_local = np.cross([0,1,0], x_local)
    y_local /= np.linalg.norm(y_local)
    z_local = np.cross(x_local, y_local)
    return np.column_stack((x_local, y_local, z_local))

def angle_between(v1, v2):
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    cos_angle = np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
    return np.degrees(np.arccos(cos_angle))

# =====================
# CÁLCULO DE POSICIONES DEFINITIVAS CON TILT COMPLETO
# =====================
crystals_positions_def = {}
tilt_contributions = {}
clover_keys = list(clovers.keys())

for idx, clover_label in enumerate(clover_keys):
    theta, phi, r_center = clovers[clover_label]

    # Centro nominal del clover
    center_nom = spherical_to_cartesian(theta, phi, r_center, depth=b)
    n_nom = center_nom / np.linalg.norm(center_nom)

    # Normal medida (end_vector)
    end_vec = end_vectors[idx]
    n_meas = end_vec / np.linalg.norm(end_vec)

    # --- Tilt: rotacional + translacional ---
    tilt_rot = angle_between(n_nom, n_meas)
    s = np.dot(center_nom, n_meas)
    center_proj = s * n_meas
    delta_vec = center_proj - center_nom  # traslación equivalente del tilt

    tilt_contributions[clover_label] = {
        'tilt_rot_deg': tilt_rot,
        'delta_vector_cm': delta_vec,
        'tilt_trans_eq_deg': np.degrees(np.arctan2(np.linalg.norm(delta_vec - np.dot(delta_vec,n_nom)*n_nom), np.linalg.norm(center_nom)))
    }

    # Base local nominal y con tilt
    M_nom = build_local_basis_face(n_nom)
    r_tilted = n_meas  # la nueva normal es la medida
    M_tilt = build_local_basis_face(r_tilted)

    # Roll sobre la normal medida
    roll_rad = np.radians(roll_angle[idx])
    R_roll = rotation_matrix(r_tilted, roll_rad)
    M_def = M_tilt @ R_roll

    # Posiciones definitivas de los cristales
    for color in original_colors:
        offset_local = crystal_offsets[color]
        pos_def = center_nom + delta_vec + M_def @ offset_local
        crystals_positions_def[(clover_label, color)] = pos_def

# =====================
# GUARDAR POSICIONES
# =====================
with open('crystal_positions_def_final.txt','w') as f:
    for key,pos in crystals_positions_def.items():
        clover_label,color = key
        f.write(f"('{clover_label}', '{color}'): {pos[0]:.6f}, {pos[1]:.6f}, {pos[2]:.6f}\n")

# =====================
# CÁLCULO DE ANGULOS ENTRE TODOS LOS PARES
# =====================
pairs_def = []
crystal_keys = list(crystals_positions_def.keys())
for i,key1 in enumerate(crystal_keys):
    for key2 in crystal_keys[i+1:]:
        v1 = crystals_positions_def[key1]
        v2 = crystals_positions_def[key2]
        angle = angle_between(v1,v2)
        pairs_def.append((key1,key2,angle))

# Guardar ángulos
with open('angles_pairs_def_final.txt','w') as f:
    for (cl1,col1),(cl2,col2),ang in pairs_def:
        f.write(f"('{cl1}', '{col1}') - ('{cl2}', '{col2}'), {ang:.4f} \n")
        f.write(f"('{cl2}', '{col2}') - ('{cl1}', '{col1}'), {ang:.4f} \n")
# =====================
# AGRUPAR ANGULOS POR BINS DE 15°
# =====================
bins = np.arange(0,195,15)
bin_counts = {f"{bins[i]}-{bins[i+1]}": [] for i in range(len(bins)-1)}

for key1,key2,angle in pairs_def:
    for i in range(len(bins)-1):
        if bins[i] <= angle < bins[i+1]:
            bin_counts[f"{bins[i]}-{bins[i+1]}"].append((key1,key2,angle))
            break

# Guardar combinaciones por bin
with open("combinaciones_por_bin_def_final.txt","w") as f:
    f.write("Bin_angular (deg) : Número de combinaciones : Pares (cristal1, cristal2, ángulo)\n\n")
    for bin_label, pairs in bin_counts.items():
        f.write(f"{bin_label} : {len(pairs)*2} combinaciones\n")  # contar duplicados
        for k1, k2, ang in pairs:  # mostrar solo pares únicos
            f.write(f"    {k1} - {k2} , {ang:.2f} deg\n")
            f.write(f"    {k2} - {k1} , {ang:.2f} deg\n")

        f.write("\n")
