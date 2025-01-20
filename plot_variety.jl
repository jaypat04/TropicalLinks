using Oscar
using Plots

# Step 1: Define a polynomial ring
R, (x,y) = polynomial_ring(QQ, ["x", "y"])

# Step 2: Define an ideal & valuation
I = ideal([x - 2 , y - 8])  # Example ideal
nu = tropical_semiring_map(QQ, 2) # p-adic valuation

# Step 3: Compute the tropical variety
trop_var = tropical_variety(I, nu; weighted_polyhedral_complex_only=true)  # Compute the tropical variety

# Step 4: Extract vertices
verts = vertices(first(trop_var))  # Get vertices of the tropical variety
verts_coords = [Tuple(v) for v in verts]  # Convert PointVector to tuples

# Step 5: Prepare data for plotting
x_coords = [Float64(c[1]) for c in verts_coords]
y_coords = [Float64(c[2]) for c in verts_coords]

# Step 6: Plot using Plots.jl
scatter(x_coords, y_coords, label="Tropical Variety Vertices", legend=:top)