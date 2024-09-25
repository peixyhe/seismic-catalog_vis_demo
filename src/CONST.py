"""
Constants and Common Functions Definition Script

This script defines constants related to a geographic region and computes conversion factors between degrees and kilometers/meters,
adjusted for the Earth's curvature at the average latitude of the region.

"""

import math

# Earth's radius in kilometers and meters
EARTH_RADIUS_KM = 6_371.393     # in kilometers
EARTH_RADIUS_M = 6_371_393.0    # in meters

# Geographic region boundaries (in degrees)
# Latitude ranges from -90 to 90 degrees
latitude_up = 33.4      # Northern boundary latitude
latitude_down = 21.8    # Southern boundary latitude

# Longitude ranges from -180 to 180 degrees
longitude_left = 97.8   # Western boundary longitude
longitude_right = 107.0 # Eastern boundary longitude

# Other constants
Z_SCALE = 10            # Scaling factor for Z-axis
MIN_MAGNITUDE = 1       # Minimum magnitude for seismic events

# Calculate the average latitude of the region
average_latitude = (latitude_up + latitude_down) / 2.0

# Calculate the conversion factor from degrees to kilometers
# This formula computes the distance corresponding to one degree difference in coordinates,
# adjusted by the scaling factor and Earth's curvature at the given latitude.
# It accounts for the fact that the length of a degree of longitude varies with latitude.

# Conversion factor from degrees to kilometers
angle_to_kilometers = (
    (2 * math.pi * EARTH_RADIUS_KM) / 360.0 / Z_SCALE
) * (
    math.sqrt(1 + math.cos(math.radians(average_latitude)) ** 2) / math.sqrt(2)
)

# Conversion factor from degrees to meters
angle_to_meters = (
    (2 * math.pi * EARTH_RADIUS_M) / 360.0 / Z_SCALE
) * (
    math.sqrt(1 + math.cos(math.radians(average_latitude)) ** 2) / math.sqrt(2)
)

print(f"One degree is equal to {angle_to_kilometers} kilometers")
print(f"Z-axis data scale: {-1.0 / angle_to_kilometers}\n")


def Get_mag2radius(mag0):
    # Convert magnitude to radius using a power function
    return round(math.pow(2.0, (mag0 - 4.56) / 1.96), 6)

def geodetic_to_ecef_m(lon_deg, lat_deg, depth_m):
    # WGS84 ellipsoid parameters
    a = 6378137.0  # Semi-major axis, in meters
    f = 1 / 298.257223563  # Flattening factor
    e2 = 2 * f - f ** 2  # Square of the first eccentricity

    # Convert angles from degrees to radians
    lon_rad = math.radians(lon_deg)
    lat_rad = math.radians(lat_deg)

    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)
    sin_lon = math.sin(lon_rad)
    cos_lon = math.cos(lon_rad)

    # Calculate the radius of curvature in the prime vertical N
    N = a / math.sqrt(1 - e2 * sin_lat ** 2)

    # Calculate ECEF coordinates (x_m, y_m, z_m)
    x_m = (N + depth_m) * cos_lat * cos_lon
    y_m = (N + depth_m) * cos_lat * sin_lon
    z_m = (N * (1 - e2) + depth_m) * sin_lat

    return x_m, y_m, z_m

def geodetic_to_ecef_km(lon_deg, lat_deg, depth_m):
    # WGS84 ellipsoid parameters
    a = 6378137.0  # Semi-major axis, in meters
    f = 1 / 298.257223563  # Flattening factor
    e2 = 2 * f - f ** 2  # Square of the first eccentricity

    # Convert angles from degrees to radians
    lon_rad = math.radians(lon_deg)
    lat_rad = math.radians(lat_deg)

    sin_lat = math.sin(lat_rad)
    cos_lat = math.cos(lat_rad)
    sin_lon = math.sin(lon_rad)
    cos_lon = math.cos(lon_rad)

    # Calculate the radius of curvature in the prime vertical N
    N = a / math.sqrt(1 - e2 * sin_lat ** 2)

    # Calculate ECEF coordinates (x_m, y_m, z_m)
    x_m = (N + depth_m) * cos_lat * cos_lon
    y_m = (N + depth_m) * cos_lat * sin_lon
    z_m = (N * (1 - e2) + depth_m) * sin_lat

    return x_m / 1000.0, y_m / 1000.0, z_m / 1000.0
