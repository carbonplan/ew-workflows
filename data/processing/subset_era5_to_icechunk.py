
import xarray as xr
from distributed import Client
from zarr.storage import ObjectStore
import icechunk
from obstore.store import GCSStore 
import zarr
from distributed import Client 
from icechunk.xarray import to_icechunk

client = Client(n_workers=64)
client 
zarr.config.set({'async.concurrency': 128})

gcs_store = GCSStore.from_url(
    'gs://gcp-public-data-arco-era5/ar/full_37-1h-0p25deg-chunk-1.zarr-v3/',skip_signature=True
)

zarr_store = ObjectStore(store= gcs_store, read_only=True)
ds = xr.open_zarr(zarr_store, consolidated=False)
varlist = ["2m_temperature",
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

ds_subset = ds[varlist]

minlat, maxlat = 24, 50
minlon, maxlon = -125, -65
mintime, maxtime = '2000', '2020'

ds_subset = ds_subset.sel(latitude=slice(maxlat, minlat), 
                                longitude=slice(360+minlon, 360+maxlon),
                                time=slice(mintime, maxtime),
                           level = 1000,
                            )

storage = icechunk.s3_storage(bucket="carbonplan-carbon-removal", prefix="era5/preprocessed_icechunk", from_env=True)
repo = icechunk.Repository.create(storage)

for var in list(ds_subset):
    print(f'writing {var}')
    ds_var = ds_subset[[var]].chunk({'time': -1, 'latitude': 12, 'longitude':12}).drop_encoding()
    
    session = repo.writable_session("main")
    to_icechunk(ds_var, session, mode='a')

    snapshot = session.commit(f"{var}")
    print(snapshot)


## Roundtrip dataset
# repo = icechunk.Repository.open(storage)
# session = repo.readonly_session("main")
# rtds = xr.open_zarr(session.store, consolidated=False)

