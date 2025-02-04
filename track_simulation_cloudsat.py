import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def plot_cloudsat_ground_track(num_orbits=10, projection_type='Orthographic', central_longitude=0, central_latitude=0):
    """
    Plots the CloudSat satellite ground track.
    
    Parameters:
    - num_orbits (int): Number of complete orbits to plot.
    - projection_type (str): Map projection type ('Orthographic', 'Mercator', etc.).
    - central_longitude (float): Center longitude for the map.
    - central_latitude (float): Center latitude for the map.
    """
    # Constants
    R_EARTH = 6371  # Radius of Earth in km
    NUM_POINTS = 10000  # More points for smoother trajectory
    EARTH_ROTATION_RATE = 360 / (24 * 60 * 60)  # Degrees per second
    
    # CloudSat Orbital Parameters
    a = 7078  # Semi-major axis in km (Earth's radius + 705 km altitude)
    e = 0.001  # Almost circular orbit
    i = np.radians(98.2)  # Inclination in radians (Sun-synchronous)
    omega = np.radians(0)  # Argument of Periapsis (set to 0)
    
    # Orbital period using Kepler's 3rd Law (T in seconds)
    T = 2 * np.pi * np.sqrt(a**3 / 398600)  # Earth's standard gravitational parameter = 398600 km^3/s^2
    
    # Time array for specified orbits
    t = np.linspace(0, num_orbits * T, NUM_POINTS)  # Time in seconds
    
    # Mean anomaly (M)
    M = (2 * np.pi / T) * t  # Linear increase over time
    
    # Eccentric anomaly (E) using Kepler's equation (approximation)
    E = M + e * np.sin(M)  # Simplified for small e
    
    # True anomaly (ν)
    nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))
    
    # Distance from Earth's center to the satellite
    r = a * (1 - e * np.cos(E))
    
    # RAAN (Ω) progression to simulate Earth rotation (degrees)
    Omega = np.radians(EARTH_ROTATION_RATE * t)
    
    # Orbital coordinates in the inertial frame
    x_orbital = r * (np.cos(nu) * np.cos(omega) - np.sin(nu) * np.sin(omega) * np.cos(i))
    y_orbital = r * (np.cos(nu) * np.sin(omega) + np.sin(nu) * np.cos(omega) * np.cos(i))
    z_orbital = r * (np.sin(nu) * np.sin(i))
    
    # Rotate coordinates to account for RAAN (Ω)
    x_eci = x_orbital * np.cos(Omega) - y_orbital * np.sin(Omega)
    y_eci = x_orbital * np.sin(Omega) + y_orbital * np.cos(Omega)
    z_eci = z_orbital  # No change in z-direction
    
    # Convert to latitude and longitude
    latitudes = np.degrees(np.arcsin(z_eci / np.sqrt(x_eci**2 + y_eci**2 + z_eci**2)))
    longitudes = np.degrees(np.arctan2(y_eci, x_eci))
    
    # Ensure longitude is within [-180, 180]
    longitudes = (longitudes + 180) % 360 - 180
    
    # Set projection
    projection_dict = {
        'Orthographic': ccrs.Orthographic(central_longitude=central_longitude, central_latitude=central_latitude),
        'Mercator': ccrs.Mercator(),
        'PlateCarree': ccrs.PlateCarree(),
        'Robinson': ccrs.Robinson(),
        'Mollweide': ccrs.Mollweide(),
        'LambertConformal': ccrs.LambertConformal(),
    }
    
    projection = projection_dict.get(projection_type, ccrs.Orthographic(central_longitude, central_latitude))
    
    # Plot using Cartopy
    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(projection=projection)
    # ax.set_global()
    ax.coastlines()
    
    # Plot the ground track
    ax.plot(longitudes, latitudes, color='#4C6EE6', linewidth=1, transform=ccrs.Geodetic(), alpha=1)
    
    plt.title(f'CloudSat Satellite Ground Track ({num_orbits/14.6} Day(s))', pad=30)
    plt.savefig(f'CloudSat_Track_days-{num_orbits/14.6}.jpg')
    plt.show()



# Example usage
# plot_cloudsat_ground_track(num_orbits=1*14.6, projection_type='Orthographic', central_longitude=75, central_latitude=35)


days = [1,2,16,30]
for i in days:
    plot_cloudsat_ground_track(num_orbits=i*14.6, projection_type='Orthographic', central_longitude=75, central_latitude=35)
