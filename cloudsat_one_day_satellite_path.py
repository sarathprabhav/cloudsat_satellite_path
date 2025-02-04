import os
from pyhdf.SD import SDC
from pyhdf.HDF import *
from pyhdf.VS import *
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


# Folder containing HDF files
data_folder = "Data"
data_folder = "/media/sarath-prabhavu/Achu/Data_cloudsat/2B_GEOPROF_LIDAR/2007/Apr"
day_name = "2007093"

# Lists to store latitude and longitude values
all_latitudes = []
all_longitudes = []

# Loop through all files in the "Data" folder
for file_name in os.listdir(data_folder):
    if file_name.endswith(".hdf") and file_name.startswith(day_name):  # Ensure it's an HDF file and starts with day_name
        file_path = os.path.join(data_folder, file_name)

        try:
            # Open the HDF file
            f = HDF(file_path, SDC.READ)
            vs = f.vstart()

            # Read latitude and longitude data
            vdata_lat = vs.attach('Latitude')
            vdata_lon = vs.attach('Longitude')

            lat = vdata_lat[:]
            lon = vdata_lon[:]

            # Append data to lists
            all_latitudes.extend(lat)
            all_longitudes.extend(lon)

            # Detach datasets and close the file
            vdata_lat.detach()
            vdata_lon.detach()
            vs.end()
            # f.end()
        
        except Exception as e:
            print(f"Error processing {file_name}: {e}")

# Now, all_latitudes and all_longitudes contain the aggregated data
print(f"Total latitude values collected: {len(all_latitudes)}")
print(f"Total longitude values collected: {len(all_longitudes)}")



# """
# Create a figure and set up an Orthographic projection
fig = plt.figure(figsize=(10, 10))
# fig.subplots_adjust(left=0.1, bottom=0, right=0.2, top=0.1)
ax = plt.axes(projection=ccrs.Orthographic(central_longitude=110, central_latitude=35))


# Ensure the entire visible hemisphere is shown
ax.set_global()

# Add features to the globe
ax.add_feature(cfeature.LAND, edgecolor='black', color='white')
ax.add_feature(cfeature.COASTLINE)
# ax.add_feature(cfeature.BORDERS, linestyle=':')

# Plot the satellite track
plt.plot(all_longitudes, all_latitudes, color='blue', linewidth=2, transform=ccrs.Geodetic(), label='Satellite Track')

# Add title and legend
plt.title("Satellite Track on Spherical Earth Map", fontsize=16, pad=30)
plt.legend()
plt.savefig("track"+day_name+".jpg")
# Show the plot
plt.show()

# """