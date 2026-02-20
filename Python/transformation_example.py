# Signed logarithmic transformation applied to all pain threshold measurements
# to normalize distributions while preserving sign information

def signed_log_transform(data):
    """Apply signed log transformation: sign(x) × log(|x| + 1)"""
    data_min = np.nanmin(data)
    s = np.sign(data + 1)  # preserve direction
    data_abs = np.absolute(data + 1)
    data_transformed = s * np.log(data_abs.astype("float"))
    return data_transformed


# Variables transformed:
# Non-sensitized: vonFrey, Hitze (heat), Kaelte (cold), Druck (pressure), Strom (electrical)
# Sensitized: vonFrey_C, Hitze_C, Kaelte_M (with capsaicin/menthol)

# Sensitization effects computed as differences:
CapsHeat = log(Hitze + 1) - log(Hitze_C + 1)  # capsaicin effect on heat
CapsvFrey = log(vonFrey + 1) - log(vonFrey_C + 1)  # capsaicin effect on von Frey
MenthCold = log(Kaelte + 1) - log(Kaelte_M + 1)  # menthol effect on cold

# Output: dfPainThresholdsAnalyzed_log_eff.csv (n=125, 11 variables)