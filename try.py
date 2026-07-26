import ee
import folium

# Authenticate and initialize with the specific project ID
ee.Authenticate()
ee.Initialize(project='ee-ishansinhagzb')

import ee
import folium

# Initialize the Earth Engine library
try:
    ee.Initialize(project='ee-ishansinhagzb')
except Exception as e:
    ee.Authenticate()
    ee.Initialize(project='ee-ishansinhagzb')

# Define the dataset and filter for a specific date range
dataset = ee.ImageCollection('NASA/GRACE/MASS_GRIDS_V04/MASCON') \
                  .filter(ee.Filter.date('2002-03-31', '2024-09-30'))

# Select the equivalent water thickness band and get the mean image for that period
equivalentWaterThickness = dataset.select('lwe_thickness').mean()

# Get the USA boundary from the LSIB dataset
usa_boundary = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017') \
    .filter(ee.Filter.eq('country_na', 'United States'))

# Clip the data to the USA boundary
equivalentWaterThickness_usa = equivalentWaterThickness.clip(usa_boundary)

# Define visualization parameters
vis_params = {
  'min': -25.0,
  'max': 25.0,
  'palette': ['001137', '01abab', 'e7eb05', '620500'],
}

# Helper function to add EE layers to folium
def add_ee_layer(self, ee_image_object, vis_params, name):
  map_id_dict = ee.Image(ee_image_object).getMapId(vis_params)
  folium.raster_layers.TileLayer(
    tiles = map_id_dict['tile_fetcher'].url_format,
    attr = 'Map Data &copy; <a href="https://earthengine.google.com/">Google Earth Engine</a>',
    name = name,
    overlay = True,
    control = True
  ).add_to(self)

# Add the method to folium.Map
folium.Map.add_ee_layer = add_ee_layer

# Create a map centered on the USA
my_map = folium.Map(location=[37.0902, -95.7129], zoom_start=4)

# Add the layer to the map
my_map.add_ee_layer(equivalentWaterThickness_usa, vis_params, 'Equivalent Water Thickness (USA)')

# Add a layer control panel
my_map.add_child(folium.LayerControl())

# Display the map
