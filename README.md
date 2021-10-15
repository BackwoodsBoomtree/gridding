# Generic Gridding Routine
One single Julia code `gridL2_Dates.jl` to grid all kinds of satellite data onto rectangular grid at arbitrary spacing (and determined spatial and temporal resolution). Input defined just via json files, which makes it generic as long as you have corner coordinates for your satellite dataset (you can grid XCO2, methane CO, etc from TROPOMI or SIF from OCO-2, OCO-3, TROPOMI, GOME-2, etc).

The main program, gridL2_Dates.jl, can compute gridded averages and save everything into a netCDF4 file that can be read via tools such as python, Julia, Panoply, etc ( 3D dataset with time/lat/lon). Most importantly, it can oversample, i.e. the gridding takes the actual footprint overlap with the final grid into account. This is done by splitting a footprint via its corner coordinates into a set of points nxn within that footprint (typically n=10, i.e. 100 points). For these, a simple "in-the-box" gridding routine will be applied, i.e. tha algorithm finds in which grid box each sub-pixel of a footprint falls. This results in fractional "samples" per grid box.

The program was written as a first step into Julia (as everything else was too slow and C++ development time just too long). So, some things are still clumsy and not really "Julian". Feel free to make things more elegant and create a pull request if you have done so. 

## Install Julia

The current version works well with Julia versions > 1.0 but I didn't test all of them. Download it for your platform from [Julia's old
releases](https://julialang.org/downloads/oldreleases/#v131_dec_30_2019).

You might have to install a couple of Julia packages, this can done like in Julia using the built-in package manager (press `]` at the Julia prompt):

```julia
julia> ]
(v1.3) pkg> add ArgParse, Statistics, Glob, JSON, Dates, Printf
```

## How to run the program

Hint for first time users, you need to use the console/terminal to use it properly. On a Mac, it is good to add an alias, e.g. 

```
alias julia /Applications/Julia-1.3.app/Contents/Resources/julia/bin/julia
```

in your .bash_profile or whatever you are using.

Once this is done, you can test out what options there are (I often do this myself):

```
$ julia ./gridL2_Dates.jl --help
usage: gridL2_Dates.jl [--Dict DICT] [-o OUTFILE] [--monthly]
                       [--compSTD] [--latMin LATMIN] [--latMax LATMAX]
                       [--lonMin LONMIN] [--lonMax LONMAX]
                       [--dLat DLAT] [--dLon DLON]
                       [--startDate STARTDATE] [--stopDate STOPDATE]
                       [--byDays BYDAYS] [--modLike] [--dateCons]
                       [--permute] [--esaVIs] [-h]

optional arguments:
  --Dict DICT           JSON dictionary file to use (default:
                        "/home/boomtree/Git/Julia/gridding/jsonFiles/tropomi_all.json")
  -o, --outFile OUTFILE
                        output filename (default TROPOMI_SIF_map.nc)
                        (default: "TROPOMI_SIF_map.nc")
  --monthly             Use time-steps in terms of months (not days)
  --compSTD             compute standard deviation within dataset
  --latMin LATMIN       Lower latitude bound (type: Float32, default:
                        -90.0)
  --latMax LATMAX       Upper latitude bound (type: Float32, default:
                        90.0)
  --lonMin LONMIN       Lower longitude bound (type: Float32, default:
                        -180.0)
  --lonMax LONMAX       Upper longitude bound (type: Float32, default:
                        180.0)
  --dLat DLAT           latitude resolution (type: Float32, default:
                        0.2)
  --dLon DLON           longitude resolution (type: Float32, default:
                        0.2)
  --startDate STARTDATE
                        Start Date (in YYYY-MM-DD) (default:
                        "2018-03-06")
  --stopDate STOPDATE   Stop Date (in YYYY-MM-DD) (default:
                        "2018-10-31")
  --byDays BYDAYS       Time steps in days (or months if --monthly is
                        set) (type: Int64, default: 8)
  --modLike             Is temporal resolution 8-day MODIS-like?
                        (default false)
  --dateCons            Conserves the --stopDate if --byDays interval
                        forces a date range that extends beyond the
                        --stopDate. (default false)
  --permute             Permute output dataset? This reorders the
                        dimensions to the conventional order of
                        time,lat,lon (z,y,x). Must have nco installed
                        in your system. (default false)
  --esaVIs              Only for the ESA TROPOMI SIF product. Grids
                        NDVI and NIRv. json file must contain keys for
                        REF_665 and REF_781.
  -h, --help            show this help message and exit
```

The most important part of this script is `--Dict` as most description about the files you want to grid are within the json file you will use here. You can also change the spatial resolution (e.g. `--dLat`), the lat/lon boundaries if you don't want the globe and the date ranges and temporal resolution you are interested in (`--startDate`, `--stopDate`, `--byDays` defines the number of days you want to aggregate together (if you set `--monthly`, this will be the step in months). If you set `--compSTD`, the standard deviation of the data within a grid box will be computed (the original variable name will be kept with an `_std` added to it and currently discarding all attributes). This is not necessarily an uncertainty estimate but the true spread of data within a grid box. 

An example on how to run this is

```
julia ./gridL2_Dates.jl --Dict ./jsonFiles/tropomi_esa_sif.json -o TROPOMI.SIF.2018-2020.0125deg.modisLike.esa.nc --startDate 2018-05-01 --stopDate 2021-08-31 --byDays 8 --dLat 0.125 --dLon 0.125 --lonMin -125 --lonMax -67 --latMin 25 --latMax 53 --modLike --dateCons --esaVIs --permute
```

### JSON files
Now to the most important part, the json structure makes the code very general (check out https://www.json.org/json-en.html).
An example file content is below:
```
{
	"basic":{
		"lat"      : "PRODUCT/latitude",
    	"lon"      : "PRODUCT/longitude",
    	"lat_bnd"  : "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds",
    	"lon_bnd"  : "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds",
    	"sza"      : "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/solar_zenith_angle",
    	"vza"      : "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/viewing_zenith_angle",
    	"time"     : "PRODUCT/time",
		"filePattern" : "TROPOSIF_L2B_YYYY-MM-DD",
		"folder"      : "/mnt/g/TROPOMI/esa/sif/v2.1/l2b"
},
	"grid":{
		"SIF_735"      : "PRODUCT/SIF_735",
		"SIF_743"      : "PRODUCT/SIF_743",
	    "SIF_Corr_735" : "PRODUCT/SIF_Corr_735",
	    "SIF_Corr_743" : "PRODUCT/SIF_Corr_743",
	    "SIF_ERROR_735": "PRODUCT/SIF_ERROR_735",
	    "SIF_ERROR_743": "PRODUCT/SIF_ERROR_743",
        "REF_665" : "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/TOA_RFL",
        "REF_680" : "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/TOA_RFL",
        "REF_712" : "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/TOA_RFL",
        "REF_741" : "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/TOA_RFL",
        "REF_755" : "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/TOA_RFL",
        "REF_773" : "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/TOA_RFL",
        "REF_781" : "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/TOA_RFL",
	    "cloud_fraction_L2" : "PRODUCT/SUPPORT_DATA/INPUT_DATA/cloud_fraction_L2"
},
    "filters":{
        "filter_lt" : {
			"PRODUCT/SUPPORT_DATA/INPUT_DATA/cloud_fraction_L2": 0.20,
            "PRODUCT/SIF_735" : 5.0,
            "PRODUCT/SIF_Corr_735" : 3.0,
            "PRODUCT/SIF_ERROR_735" : 0.3
		}
},
	"glob_attribs":{
		"title"       : "TROPOMI SIF and Reflectances as Derived from ESA L2B Product",
		"institution" : "University of Oklahoma, College of Atmospheric and Geographic Sciences",
		"author"      : "Russell Doughty, PhD",
		"source"      : "This gridded data was generated from the ESA TROPOMI L2B SIF product available at https://s5p-troposif.noveltis.fr/data-access/.",
		"references"  : "Guanter, L., Bacour, C., Schneider, A., Aben, I., van Kempen, T.A., Maignan, F., Retscher, C., Köhler, P., Frankenberg, C., Joiner, J. and Zhang, Y., 2021. The TROPOSIF global sun-induced fluorescence dataset from the Sentinel-5P TROPOMI mission. Earth System Science Data Discussions, pp.1-27.",
		"comment"     : "Gridding was done with the oversampling method, originally coded by Dr. Christian Frankenberg (https://github.com/cfranken/gridding). Data was permuted using the NCO package.",
        "filters"     : "Input soundings were filtered for cloud fracton < 0.2, SIF_735 < 5.0, SIF_Corr_735 < 3.0, and SIF_ERROR_735 < 0.3."
}
}
```

The structure will be used to define dictionaries in Julia, with a key (the left side, pointing to internal variables in the code) and a value (this is actually the path to the respetive dataset within the files you want to grid). What is needed for sure is the key/value pair for `lat_bnd` and `lon_bnd` as these corner coordinates in the files are used to perform proper gridding and over-sampling (I think time is not even used right now). The other part is in the `grid` group, eveything you will add here will be gridded, i.e. all variables on the left side of the key/value pairs (and saved how you call it on the left). 

Last but not least, we need a `filePattern`, which defines where to look for the data (`YYYY`, `MM` and `DD` are keywords, which will be internally used to find the right matching years, months and days). Wildcards, such as `?` and `*` are not necessary as patial matches are returned automatically. Then you need to provide the main folder where the data is located (the full path is `folder/filePattern`). All files, including those in subdirectories, will be scanned using the `filePattern`.

### Filter criteria:
You can add filter criteria as well, as we sometimes want to make sure that quality filters are applied, angles are within a specific range, and so forth. Within the json file, this can be done by adding groups called `filter_eq`, `filter_gt`, or `filter_lt` in the `filters` dictionary. These test for equalities, greater than, or lower than. Here is an example:
```
  "filters":{
      "filter_lt" : {
          "PRODUCT/SUPPORT_DATA/INPUT_DATA/cloud_fraction_L2": 0.20,
          "PRODUCT/SIF_735" : 5.0,
          "PRODUCT/SIF_Corr_735" : 3.0,
          "PRODUCT/SIF_ERROR_735" : 0.3
        }
    }
```
Here, we want to make sure that the `cloud_fraction` is less than than 0.20 and that the SIF values are within a reasonable range (be careful doing things like this, just added for demonstration).

## Code of Conduct:
Please feel free to use this tool but make sure that you help the community if you find bugs, improve it etc. Any modifications that are useful should be made publicly available, you can fork and create a pull request. Also, let us know if you find bugs. On top of that, please acknowledge the tool if you use it in publications.

## MODIS files
We use these as well but I have no time to document all of them right now. Contact us if you want to use them (cfranken"at"caltech.edu).
