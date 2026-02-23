





import xarray as xr
from zarr.storage import ObjectStore
import icechunk
from obstore.store import GCSStore 
import zarr
from icechunk.xarray import to_icechunk
import coiled 
import flox
from scipy import stats

zarr.config.set({'async.concurrency': 128})


cluster = coiled.Cluster(n_workers=[30,50], name='era5',worker_vm_types = ['c8g.2xlarge','m8g.2xlarge'], scheduler_vm_types='c8g.xlarge',spot_policy='spot_with_fallback', region='us-west-2',tags={'Project':'ERW'})
client = cluster.get_client()
print(client)


varlist = [
    "2m_temperature",
        "skin_reservoir_content",
        "volumetric_soil_water_layer_1",
        "volumetric_soil_water_layer_2",
        "volumetric_soil_water_layer_3",
        "volumetric_soil_water_layer_4",
        "soil_temperature_level_1",
        "soil_temperature_level_2",
        "soil_temperature_level_3",
        "soil_temperature_level_4",
        "potential_evaporation",
        "runoff",
        "surface_runoff",
        "sub_surface_runoff",
        "evaporation",
        "total_precipitation",
        "geopotential",
        "land_sea_mask",
        "soil_type"
        ]
        
mintime, maxtime = '2000', '2020'

gcs_store = GCSStore.from_url(
    'gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3/',skip_signature=True
)
zarr_store = ObjectStore(store= gcs_store, read_only=True)
ds = xr.open_zarr(zarr_store, consolidated=False, chunks={'time':30,'latitude':721,'longitude':1440}).drop_encoding()
ds = ds[varlist]
ds = ds.sel(time=slice(mintime, maxtime),level = 1000).drop_vars('level')
ds.coords['longitude'] = (ds.coords['longitude'] + 180) % 360 - 180
ds_subset = ds.sortby(ds['longitude'])
ds_subset.attrs = {}



vars_mean = [v for v in ds_subset.data_vars if v != 'soil_type']
# coarsen roughly 1/4 degree to 1 degree

ds_coarse = ds_subset[vars_mean].coarsen(latitude=4, longitude=4, boundary='trim').mean()

# for soil_type, we're doing mode instead of mean for our coarsening. Using scipy + xarray apply_ufunc to keep it lazy. 
soil_constructed = ds_subset['soil_type'].isel(time=0).coarsen(
    latitude=4, longitude=4, boundary='trim'
).construct(
    latitude=('latitude', 'lat_bins'), 
    longitude=('longitude', 'lon_bins')
)

def mode_func(x):
    flat = x.reshape(*x.shape[:-2], -1)
    return stats.mode(flat, axis=-1, keepdims=False).mode

soil_coarse = xr.apply_ufunc(
    mode_func,
    soil_constructed,
    input_core_dims=[['lat_bins', 'lon_bins']],
    vectorize=False,
    dask='parallelized',
    output_dtypes=[ds_subset['soil_type'].dtype]
)


ds_coarse['soil_type'] = soil_coarse


# setup chunk/shard encoding
encoding = {}
  
for var in list(ds):
    # soil_type is reduced to lat/lon, so we don't need time chunking
    if var =='soil_type':
        encoding[var] =     {
        "chunks": (3,3),     
        "shards": (21,21),      
    }
    else:
        encoding[var] =     {
            "chunks": (184104, 3,3),     
            "shards": (184104, 21,21),      
        }


# chunk to shard size
rds = ds_coarse.chunk({'time': -1, 'latitude': 21, 'longitude':21}) 
rds['soil_type'] = rds['soil_type'].chunk({'latitude': -1, 'longitude':-1})

storage = icechunk.s3_storage(bucket="carbonplan-carbon-removal", prefix="era5/preprocessed_icechunk", from_env=True)
repo = icechunk.Repository.open_or_create(storage)

session = repo.writable_session("main")
to_icechunk(rds, session, encoding=encoding,mode='w')

snapshot = session.commit("era5_coarsen")

client.shutdown()



## roundtrip the data
# storage = icechunk.s3_storage(bucket="carbonplan-carbon-removal", prefix="era5/preprocessed_icechunk", from_env=True)
# repo = icechunk.Repository.open(storage)

# session = repo.readonly_session("main")
# rtds = xr.open_dataset(session.store, engine='zarr', chunks="auto")
# rtds