# Examples folder inputs

This folder contains example inputs for running **uraster** workflows.  
Each example includes a source mesh, one or more raster files, and configuration details that defines how raster values are mapped or aggregated onto the mesh.

## Format used for each example

For consistency, every example entry uses the same structure.

* Source mesh, what mesh is used and its approximate resolution  
* Input raster, what raster is used, key raster metadata, and data sources
* Configuration, what key settings it must contain  

## Input mesh file source

### DGGS meshes
All DGGS mesh files used in the examples are generated using the [**DGGAL**](https://github.com/ecere/dggal/tree/eC-core) library.  

### TIN meshes
Place holder.

### MPAS meshes
Place holder.

---

## Example 1

### Source mesh
* Mesh type: DGGS
* Mesh file or files `rhealpix_global_res3.geojson`
* Cell size: 116744 km²

### Input raster
* Raster file or files `EDGAR_CH4_emission_global_2015.tiff`
* Tiling: single
* Data type: continuous
* Spatial resolution: 0.1°
* CRS: EPSG 4326
* Spatial coverage: global
* Notes: EDGAR v8.0 Methane Emissions inventory data from manure management sector for the year 2015, unit: ton
* Data source:  
  Data can be downloaded from [here](https://data.jrc.ec.europa.eu/dataset/b54d8149-2864-4fb9-96b9-5fd3a020c224)
* Reference:  
  Crippa, Monica; Guizzardi, Diego; Pagani, Federico; Banja, Manjola; Muntean, Marilena; Schaaf, Edwin; Becker, William; Monforti Ferrario, Fabio; Grassi, Giacomo; Rossi, Simone; Brandao De Melo, Joana; Jacome Felix Oom, Duarte; Branco, Alfredo; San-Miguel, Jesus; Vignati, Elisabetta (2023): EDGAR v8.0 Greenhouse Gas Emissions. European Commission, Joint Research Centre (JRC) [Dataset] doi: 10.2905/b54d8149-2864-4fb9-96b9-5fd3a020c224 

### Configuration
* Required fields
  * Mesh path and raster path
  * Output path and format

---

## Example 2

### Source mesh
* Mesh type: DGGS
* Mesh file or files `isea3h_bbox_res15.geojson`
* Cell size: 3.55 km²

### Input raster
* Raster file or files `WSF2015_v2_-80_42.tif`
* Tiling: single
* Data type: discrete
* Spatial resolution: 0.32 arcsec
* CRS: EPSG 4326
* Spatial coverage: one tile covering Toronto, Canada and Buffalo, US
* Notes: World Settlement Footprint 2015 data, binary built up layer where pixel value 255 indicates built up areas and 0 indicates non built up areas
* Data source:  
  Data can be downloaded from [here](https://download.geoservice.dlr.de/WSF2015/)
* Reference:  
    Marconcini, M., Metz-Marconcini, A., Üreyen, S., Palacios-Lopez, D., Hanke, W., Bachofer, F., Zeidler, J., Esch, T., Gorelick, N., Kakarla, A. and Paganini, M., 2020. Outlining where humans live, the World Settlement Footprint 2015. Scientific Data, 7(1), 242.

### Configuration
* Required fields
  * Mesh path and raster path
  * Output path and format

---

## Example 3

### Source mesh
* Mesh type: DGGS
* Mesh file or files `rhealpix_China_res3.geojson`, `rhealpix_China_res6.geojson`
* Cell size: res3 116744 km², res6 160 km²

### Input raster
* Raster file or files `China_CH4_emission_2020.tif`
* Tiling: single
* Data type: continuous
* Spatial resolution: 10 km
* CRS: EPSG 4024
* Spatial coverage: China
* Notes: Methane emissions from coal mining in China for the year 2020, unit: Mg km⁻² a⁻¹ 
* Data source:  
  Data can be requested via [link](https://forms.gle/NGMXUTfMumMFkMZPA)
* Reference:  
    Sheng, J., Song, S., Zhang, Y., Prinn, R.G. and Janssens-Maenhout, G., 2019. Bottom-up estimates of coal mine methane emissions in China: A gridded inventory, emission factors, and trends. Environmental Science & Technology Letters, 6(8), pp.473-478.

### Configuration
* Required fields
  * Mesh path and raster path
  * Output path and format

---

## Example 4

### Source mesh
* Mesh type: DGGS
* Mesh file or files `isea7h_bbox_res5.geojson`
* Cell size: 3034.84 km²

### Input raster
* Raster file or files  
  `Hansen_GFC-2019-v1.7_treecover2000_40N_090W.tif`,  
  `Hansen_GFC-2019-v1.7_treecover2000_40N_100W.tif`,  
  `Hansen_GFC-2019-v1.7_treecover2000_50N_090W.tif`,  
  `Hansen_GFC-2019-v1.7_treecover2000_50N_100W.tif`
* Tiling: multi tile
* Data type: continuous
* Spatial resolution: 0.00025°
* CRS: EPSG 4326
* Spatial coverage: four tiles covering the inland eastern United States
* Notes: Tree cover in the year 2000, defined as canopy closure for all vegetation taller than 5 m in height. Encoded as a percentage per output grid cell, with values ranging from 0 to 100
* Data source:  
  Data can be downloaded from [here](https://earthenginepartners.appspot.com/science-2013-global-forest/download_v1.7.html)
* Reference:  
    Hansen, M.C., Potapov, P.V., Moore, R., Hancher, M., Turubanova, S.A., Tyukavina, A., Thau, D., Stehman, S.V., Goetz, S.J., Loveland, T.R. and Kommareddy, A., 2013. High-resolution global maps of 21st-century forest cover change. science, 342(6160), pp.850-853.

### Configuration
* Required fields
  * Mesh path and raster path
  * Output path and format

---

## Example 5

### Source mesh
* Mesh type: DGGS
* Mesh file or files `isea7h_bbox_res8.geojson`
* Cell size: 1.26 km²

### Input raster
* Raster file or files  
  `ArcticDEM_32_34_1_1_2m_v4.1_dem.tif`,  
  `ArcticDEM_32_34_1_2_2m_v4.1_dem.tif`,  
  `ArcticDEM_32_35_1_1_2m_v4.1_dem.tif`,  
  `N79W078_FABDEM_V1-0.tif`,  
  `N79W079_FABDEM_V1-0.tif`,  
  `N79W080_FABDEM_V1-0.tif`,  
  `N79W081_FABDEM_V1-0.tif`,  
  `N79W082_FABDEM_V1-0.tif`,  
  `N79W083_FABDEM_V1-0.tif`
* Tiling: multi tile
* Data type: continuous
* Spatial resolution:  
  ArcticDEM 2 m  
  FABDEM 1 arcsec
* CRS:  
  ArcticDEM EPSG 3413  
  FABDEM EPSG 4326
* Spatial coverage: three ArcticDEM tiles and six FABDEM tiles covering high latitude areas in Nunavut
* Notes:  
  ArcticDEM is a high resolution, high quality digital surface model of the Arctic. The DSM can be treated as a DEM in polar regions where vegetation and built structures are minimal. 
  FABDEM is a global elevation product derived from Copernicus DEM that removes building and tree height biases.
* Data source:  
  ArcticDEM data can be downloaded from [here](https://www.arcgis.com/apps/webappviewer/index.html?id=aff5fa8f5d5548c6bff44cc8be385f61)  
  FABDEM data can be downloaded from [here](https://data.bris.ac.uk/data/dataset/25wfy0f9ukoge2gs7a5mqpq2j7)

* Reference:  
    Hawker, L., Uhe, P., Paulo, L., Sosa, J., Savage, J., Sampson, C. and Neal, J., 2022. A 30 m global map of elevation with forests and buildings removed. Environmental Research Letters, 17(2), 024016.
    Porter, Claire; Howat, Ian; Noh, Myoung-Jon; Husby, Erik; Khuvis, Samuel; Danish, Evan; Tomko, Karen; Gardiner, Judith; Negrete, Adelaide; Yadav, Bidhyananda; Klassen, James; Kelleher, Cole; Cloutier, Michael; Bakker, Jesse; Enos, Jeremy; Arnold, Galen; Bauer, Greg; Morin, Paul, (2023): ArcticDEM - Mosaics, Version 4.1, Harvard Dataverse, V1 [Dataset] doi:10.7910/DVN/3VDC4W

### Configuration
* Required fields
  * Mesh path and raster path
  * Output path and format

---

## Example 6

### Source mesh
* Mesh type: DGGS
* Mesh file or files `isea3h_bbox_res12.geojson`
* Cell size: 95.98 km²

### Input raster
* Raster file or files  
  `11T_20240101-20241231.tif`,  
  `11U_20240101-20241231.tif`,  
  `12T_20240101-20241231.tif`,  
  `12U_20240101-20241231.tif`
* Tiling: multi tile
* Data type: discrete
* Spatial resolution: 10 m
* CRS: EPSG 32612, EPSG 32611
* Spatial coverage: four tiles covering part of the Rocky Mountains and surrounding areas
* Notes: Esri land use and land cover product derived from ESA Sentinel 2 imagery for the year 2024
* Data source:  
  Data can be downloaded from [here](https://livingatlas.arcgis.com/landcoverexplorer/#mapCenter=-117.27100%2C33.94200%2C10.00&mode=step&timeExtent=2017%2C2022&year=2022)
* Reference:  
    Karra, K., Kontgis, C., Statman-Weil, Z., Mazzariello, J.C., Mathis, M. and Brumby, S.P., 2021, July. Global land use/land cover with Sentinel 2 and deep learning. In 2021 IEEE international geoscience and remote sensing symposium IGARSS (pp. 4704-4707). IEEE.

### Configuration
* Required fields
  * Mesh path and raster path
  * Output path and format

---

## Example 7

### Source mesh
* Mesh type: 
* Mesh file or files 
* Cell size: 

### Input raster
* Raster file or files 
* Tiling: 
* Data type: 
* Spatial resolution: 
* CRS: 
* Spatial coverage: 
* Notes: 
* Data source:  
  Data can be downloaded from
* Reference:  

### Configuration
* Required fields
  * Mesh path and raster path
  * Output path and format
---

## Example 8

### Source mesh
* Mesh type: 
* Mesh file or files 
* Cell size: 

### Input raster
* Raster file or files 
* Tiling: 
* Data type: 
* Spatial resolution: 
* CRS: 
* Spatial coverage: 
* Notes: 
* Data source:  
  Data can be downloaded from
* Reference:  

### Configuration
* Required fields
  * Mesh path and raster path
  * Output path and format

---

## Example 9

### Source mesh
* Mesh type: TIN
* Mesh file or files 
* Cell size: 

### Input raster
* Raster file or files 
* Tiling: 
* Data type: 
* Spatial resolution: 
* CRS: 
* Spatial coverage: 
* Notes: 
* Data source:  
  Data can be downloaded from
* Reference:  

### Configuration
* Required fields
  * Mesh path and raster path
  * Output path and format

---

## Example 10

### Source mesh
* Mesh type: MPAS
* Mesh file or files 
* Cell size: 

### Input raster
* Raster file or files 
* Tiling: 
* Data type: 
* Spatial resolution: 
* CRS: 
* Spatial coverage: 
* Notes: 
* Data source:  
  Data can be downloaded from
* Reference:  

### Configuration
* Required fields
  * Mesh path and raster path
  * Output path and format

---

## Example 11

### Source mesh
* Mesh type: 
* Mesh file or files 
* Cell size: 

### Input raster
* Raster file or files 
* Tiling: 
* Data type: 
* Spatial resolution: 
* CRS: 
* Spatial coverage: 
* Notes: 
* Data source:  
  Data can be downloaded from
* Reference:  

### Configuration
* Required fields
  * Mesh path and raster path
  * Output path and format

---

## Example 12

### Source mesh
* Mesh type: 
* Mesh file or files 
* Cell size: 

### Input raster
* Raster file or files 
* Tiling: 
* Data type: 
* Spatial resolution: 
* CRS: 
* Spatial coverage: 
* Notes: 
* Data source:  
  Data can be downloaded from
* Reference:  

### Configuration
* Required fields
  * Mesh path and raster path
  * Output path and format
