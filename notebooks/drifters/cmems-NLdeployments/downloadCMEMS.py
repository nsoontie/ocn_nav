# Script to download CMEMS data
# FYI: This script identifies files on the CMEMS server to download.
# If some of those files have already been downloaded, they will be downloaded again and overwritten.
# Must be run with the copernicusmarine toolbox installed.
# Also, a credential files is needed (or credential environment variables)
# To install copernicuamaarine toolboc and set up redenctials see:
# https://help.marine.copernicus.eu/en/articles/7970514-copernicus-marine-toolbox-installation

import yaml
import copernicusmarine

dataset_id = 'cmems_obs-ins_glo_phy-cur_nrt_drifter_irr'

drifterID_file = 'NLDrifterIDs.yaml'
d = 'NLSummerAZMP2024' # The deployment to download. Should be documented in drifterID_file.

with open(drifterID_file) as f:
    deployments = yaml.safe_load(f)

if d != 'all':
    try:
        IDs = deployments[d]
    except KeyError:
        print(f'Deployment name {d} not available in {drifterID_file}. '\
               'No drifter IDs found.')
        IDs = []
else:
    IDs = []
    for d in deployments:
        IDs.extend(deployments[d])

nodata_IDs = []
for ID in IDs:
    print(f'Processing ID {ID}')
    # Plan is to start by downloading the monthly (last 5 years)
    # If empty, download latest (last 30 days)
    # If empty, download history (longest recrod)
    output =  copernicusmarine.get(
        dataset_id=dataset_id,
        show_outputnames=True,
        filter=f'*{ID}*.nc',
        force_download=True,
        overwrite_output_data=True,
        dataset_part='monthly'
    )
    if not output:
        output =  copernicusmarine.get(
            dataset_id=dataset_id,
            show_outputnames=True,
            filter=f'*{ID}*.nc',
            force_download=True,
            overwrite_output_data=True,
            dataset_part='latest'
        )
    if not output:
        output =  copernicusmarine.get(
            dataset_id=dataset_id,
            show_outputnames=True,
            filter=f'*{ID}*.nc',
            force_download=True,
            overwrite_output_data=True,
            dataset_part='history'
        )
    if not output:
        nodata_IDs.append(ID)
print(f'Finished processings IDs: {IDs}')
print(f'The following IDs contained no data: {nodata_IDs}')
