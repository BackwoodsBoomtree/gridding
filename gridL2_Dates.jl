#!/home/cfranken//julia

# Argument Parser
using ArgParse
using Base, Dates, Printf
# NetCDF tools for reading and writing
using NCDatasets
# Basic statistics
using Statistics
# File search and completion
using Glob
# JSON files
using JSON
# Parallel computing
#using Distributed, SharedArrays
# Profiler
#using Profile

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--Dict"
            help     = "JSON dictionary file to use"
            arg_type = String
            default  = "/home/boomtree/Git/Julia/gridding/jsonFiles/tropomi_all.json"
        "--outFile", "-o"
            help     = "output filename (default TROPOMI_SIF_map.nc)"
            arg_type = String
            default  = "TROPOMI_SIF_map.nc"
        "--monthly"
            help   = "Use time-steps in terms of months (not days)"
            action = :store_true
        "--compSTD"
            help   = "compute standard deviation within dataset"
            action = :store_true
        "--latMin"
            help     = "Lower latitude bound"
            arg_type = Float32
            default  = -90.0f0
        "--latMax"
            help     = "Upper latitude bound"
            arg_type = Float32
            default  = 90.0f0
        "--lonMin"
            help     = "Lower longitude bound"
            arg_type = Float32
            default  = -180.0f0
        "--lonMax"
            help     = "Upper longitude bound"
            arg_type = Float32
            default  = 180.0f0
        "--dLat"
            help     = "latitude resolution"
            arg_type = Float32
            default  = 0.20f0
        "--dLon"
            help     = "longitude resolution"
            arg_type = Float32
            default  = 0.20f0
        "--startDate"
                help     = "Start Date (in YYYY-MM-DD)"
                arg_type = String
                default  = "2018-03-06"
        "--stopDate"
                help     = "Stop Date (in YYYY-MM-DD)"
                arg_type = String
                default  = "2018-10-31"
        "--byDays"
                help     = "Time steps in days (or months if --monthly is set)"
                arg_type = Int64
                default  = 8
        "--modLike"
                help     = "Is temporal resolution 8-day MODIS-like? (default false)"
                action = :store_true
        "--dateCons"
                help     = "Conserves the --stopDate if --byDays interval forces a date range that extends beyond the --stopDate. (default false)"
                action = :store_true
        "--permute"
                help     = "Permute output dataset? This reorders the dimensions to the conventional order of time,lat,lon (z,y,x). Must have nco installed in your system. (default false)"
                action = :store_true
    end
    return parse_args(s)
end

function getFilter(name, jsonDict)
    ff = try
        jsonDict[name]
    catch
        Dict()
    end
    return ff
end

# This splits up the entire region into one grid set
function divLine!(lat1,lon1,lat2,lon2,n, points,j )
    dLat = (lat2-lat1)/(2*n)
    dLon = (lon2-lon1)/(2*n)
    startLat = lat1+dLat
    startLon = lon1+dLon
    @inbounds for i in 1:n
        #println(startLat+2*(i-1)*dLat, " ", startLon+2*(i-1)*dLon)
        #weights[(iLat-minLat+1), (iLon-minLon+1)]+=1
        points[j,i,1] = startLat+2*(i-1)*dLat
        points[j,i,2] = startLon+2*(i-1)*dLon
    end
end

# For first run (2 baselines)
function divLine2!(lat1,lon1,lat2,lon2,n, lats, lons)
    dLat = (lat2-lat1)/(2*n)
    dLon = (lon2-lon1)/(2*n)
    startLat = lat1+dLat
    startLon = lon1+dLon
    @inbounds for i in 1:n
        lats[i] = startLat+2*(i-1)*dLat
        lons[i] = startLon+2*(i-1)*dLon
    end
end

# Divide each polygon into multiple points
function getPoints!(points, vert_lat, vert_lon, n,lats_0, lons_0,lats_1, lons_1 )
    # Get reference points for two lines at the extremes:
    # println(vert_lat)
    divLine2!(vert_lat[1],vert_lon[1],vert_lat[2],vert_lon[2],n,lats_0, lons_0)
    divLine2!(vert_lat[4],vert_lon[4],vert_lat[3],vert_lon[3],n,lats_1, lons_1)
    for i in 1:n
         divLine!(lats_0[i], lons_0[i] ,lats_1[i], lons_1[i],n, points,i )
    end
end

# Still need to make sure the corners are read in properly!
function getNC_var(fin, path, DD::Bool)
    try
        loc = split(path ,r"/")
        #println(loc)
        if length(loc)==1
            return fin[path].var[:]
        elseif length(loc)>1
            gr = []
            for i in 1:length(loc)-1
                if i==1
                    gr = fin.group[loc[i]]
                else
                    gr = gr.group[loc[i]]
                end
            end

            #println(loc[end])
            si = size(gr[loc[end]])
            # DD means there is a 2nd index for footprint bounds of dimension 4!
            if DD
                if si[1]==4
                    return reshape(gr[loc[end]].var[:],4,prod(si[2:end]))'
                elseif si[end]==4
                    return reshape(gr[loc[end]].var[:],prod(si[1:end-1]),4)
                end
            else
                return reshape(gr[loc[end]].var[:],prod(si))
            end
        end
    catch e
        # @show e
        # println("Error in getNC_var ", path)
        return 0.0
    end

end

function getNC_attrib(fin, path, attri)
    loc = split(path ,r"/")
    #println(loc)
    if length(loc)==1
        return fin[path].attrib[attri]
    elseif length(loc)>1
        gr = []
        for i in 1:length(loc)-1
            if i==1
                gr = fin.group[loc[i]]
            else
                gr = gr.group[loc[i]]
            end
        end
        return gr[loc[end]].attrib[attri]
    end
end

function datesBy(startDate, endDate, byDays)
    startDate = DateTime(startDate)
    endDate   = DateTime(endDate)
    byDays    = Dates.Day(byDays)
    # Start dates for the year
    firstList = collect(startDate:byDays:endDate)
    # Last date for each range
    lastList  = firstList + byDays - Dates.Day(1)
    dateList  = hcat(firstList, lastList)
    return(dateList)
end

function modDates(startDate, endDate)
    startDate = DateTime(startDate)
    endDate   = DateTime(endDate)
    nyear     = Dates.Year(endDate).value - Dates.Year(startDate).value + 1
    yearList  = [Dates.Year(startDate).value : 1 : Dates.Year(endDate).value;]
    if nyear == 1
        # Start dates for the year
        firstList = collect(startDate:Dates.Day(8):endDate)
        # Last date for each range
        lastList  = firstList + Dates.Day(7)
        dateList  = hcat(firstList, lastList)
    elseif nyear == 2     
        # Start dates from the first year
        firstList = collect(startDate:Dates.Day(8):DateTime(string(Dates.Year(startDate).value, "-12-31")))
        # Append start dates from the second year
        firstList = [firstList;collect(DateTime(string(Dates.Year(endDate).value, "-01-01")):Dates.Day(8):endDate)]
        lastList  = firstList + Dates.Day(7)
        dateList  = hcat(firstList, lastList)
    elseif nyear > 2
        # Start dates from the first year
        firstList = collect(startDate:Dates.Day(8):DateTime(string(Dates.Year(startDate).value, "-12-31")))
        # Append start dates from middle years
        for year in yearList[1] + 1 : yearList[end] - 1
            firstList = [firstList;collect(DateTime(string(year, "-01-01")):Dates.Day(8):DateTime(string(year, "-12-31")))]
        end
        # Append start dates from final year
        firstList = [firstList;collect(DateTime(string(Dates.Year(endDate).value, "-01-01")):Dates.Day(8):endDate)]
        lastList  = firstList + Dates.Day(7)
        dateList  = hcat(firstList, lastList)
    end
    # Correct last dates to be end of year
    for year in yearList
        dateList = replace(dateList, DateTime(string(year, "-01-02")) => DateTime(string(year - 1, "-12-31")))
        dateList = replace(dateList, DateTime(string(year, "-01-03")) => DateTime(string(year - 1, "-12-31")))
    end
    return(dateList)
end

# Test for how fast we can compute point in distributed mode (not yet used)
#function compPoints!(lat,lon,n)
#    dim = size(lat)
#    lats_0 = zeros(n)
#    lons_0 = zeros(n)
#    lats_1 = zeros(n)
#    lons_1 = zeros(n)
#    a = SharedArray{Float32}((dim[1],n,n,2))
#    @distributed for i in 1:dim[1]
#        po = getPoints!(a[i,:,:,:],lat[i,:],lon[i,:],n,lats_0, lons_0,lats_1, lons_1 )
#    end
#   return a
#end

function favg_all!(arr,std_arr, weight_arr, compSTD, lat,lon,inp,s,s2,n, latMin, latMax, lonMin, lonMax, nLat, nLon, points)

    dim = size(arr)
    #println(dim)
    # Predefine some arrays to reduce allocations
    ix = zeros(Int32,n^2)
    iy = zeros(Int32,n^2)
    lats_0  = zeros(n)
    lons_0  = zeros(n)
    lats_1  = zeros(n)
    lons_1  = zeros(n)
    iLon    = floor.(Int32,lon)
    iLat    = floor.(Int32,lat)
    minLat  = minimum(Int32,floor.(lat), dims=2)
    maxLat  = maximum(Int32,floor.(lat), dims=2)
    minLon  = minimum(Int32,floor.(lon), dims=2)
    maxLon  = maximum(Int32,floor.(lon), dims=2)
    distLon = maxLon-minLon
    # How many individual grid cells might actually be there:
    dimLat = maxLat-minLat
    dimLon = maxLon-minLon
    fac    = Float32(1/n^2)

    @inbounds for i in 1:s
        #println(i, " ", dimLat[i], " ", dimLon[i])
        # Take it easy if all corners already fall into one grid box:
        if (dimLat[i]==0) & (dimLon[i]==0)
            
            weight_arr[iLon[i,1],iLat[i,1]] += 1
            for z in 1:s2
                mean_old = arr[iLon[i,1],iLat[i,1],z];
                arr[iLon[i,1],iLat[i,1],z] = mean_old .+ 1/weight_arr[iLon[i,1],iLat[i,1]] .* (inp[i,z]-mean_old);
                if compSTD
                    std_arr[iLon[i,1],iLat[i,1],z] += (inp[i,z]-mean_old) .* (inp[i,z]-arr[iLon[i,1],iLat[i,1],z])
                end
            end
        # if not, compute appropriate weights
        elseif (distLon[i])<n
            getPoints!(points,lat[i,:],lon[i,:],n,lats_0, lons_0,lats_1, lons_1 )

            ix[:] = floor.(Int32,points[:,:,1][:])
            iy[:] = floor.(Int32,points[:,:,2][:])
            
            @inbounds for j in eachindex(ix)
                weight_arr[iy[j],ix[j]] += fac;
                for z in 1:s2
                    mean_old = arr[iy[j],ix[j],z]
                    arr[iy[j],ix[j],z] = mean_old .+ fac/weight_arr[iy[j],ix[j]] .* (inp[i,z]-mean_old);
                    if compSTD
                        std_arr[iy[j],ix[j],z] +=  fac .* (inp[i,z]-mean_old) .* (inp[i,z]-arr[iy[j],ix[j],z])
                    end
                end
            end
        end
    end
end

function main()

    startTime = DateTime(now())

    #addprocs()
    # Parse command line arguments
    ar = parse_commandline()

    # Find files to be processed
    startDate = DateTime(ar["startDate"])
    stopDate  = DateTime(ar["stopDate"])
    
    if ar["monthly"]
        byDay  = Dates.Month(ar["byDays"])
    else
        byDay  = Dates.Day(ar["byDays"])
    end

    println("")
    println("Processing date range is:")
    println(Date(startDate), " to ", Date(stopDate))
    println("")

    # Build 2d array of beginning and end dates for entire range
    # If modLike but start date is bad then exit,
    # else get date ranges and count cT (timesteps)
    if ar["modLike"]
        # Exit if start date is bad
        modstartYear = Dates.Year(startDate).value
        modDateList  = collect(Date(string(modstartYear, "-01-01")):Dates.Day(8):Date(string(modstartYear, "-12-31")))
        if Date(startDate) ∉ modDateList
            println("ERROR: Input --startDate does not fall on an 8-day MODIS start date.")
            exit()
        else
            dateChunks = modDates(ar["startDate"], ar["stopDate"])
            cT = size(dateChunks, 1)
            println("Temporal resolution is modis-like. Number of time steps is: ", cT)
        end
    else
        dateChunks = datesBy(startDate, stopDate, byDay)
        cT         = size(dateChunks, 1)
        println("Temporal resolution is not modis-like. Number of time steps is: ", cT)
    end
    # Change last date according to dateCons flag
    if ar["dateCons"]
        dateChunks[end] = stopDate
        println("Warning: --dateCons true. If date range is not divisible by --byDays, then excluding data beyond the --stopDate.")
        println("")
    else
        println("Warning: --dateCons false. If date range is not divisible by --byDays, then including data beyond the --stopDate.")
        println("")
    end

    # Just lazy (too cumbersome in code as often used variables here)
    latMax = ar["latMax"]
    latMin = ar["latMin"]
    lonMax = ar["lonMax"]
    lonMin = ar["lonMin"]
    dLat   = ar["dLat"]
    dLon   = ar["dLon"]
    eps    = dLat/100

    # Define spatial grid:
    lat = collect(latMin+dLat/2.:dLat:latMax-dLat/2.0+eps)
    lon = collect(lonMin+dLon/2.:dLon:lonMax-dLon/2.0+eps)
    println("Output file dimension (time/lat/lon):")
    println(cT, "/", length(lat),"/", length(lon))
    println("")
    
    # Create output file:
    dsOut = Dataset(ar["outFile"],"c")
    defDim(dsOut,"time", cT)
    defDim(dsOut,"lat",length(lat))
    defDim(dsOut,"lon",length(lon))

    dsTime   = defVar(dsOut,"time",Float32,("time",),attrib = ["units" => "days since 1970-01-01","long_name" => "Time (UTC), start of interval"])
    dsLat    = defVar(dsOut,"lat",Float32,("lat",), attrib = ["units" => "degrees_north","long_name" => "Latitude"])
    dsLon    = defVar(dsOut,"lon",Float32,("lon",), attrib = ["units" => "degrees_east","long_name" => "Longitude"])
    dsLat[:] = lat
    dsLon[:] = lon

    # Define a global attribute
    dsOut.attrib["title"] = "TROPOMI Gridded Data"

    # Define gridded variables:
    n   = zeros(Float32,(length(lat),length(lon)))
    SIF = zeros(Float32,(length(lat),length(lon)))
    # Parse JSON files as dictionary
    jsonDict = JSON.parsefile(ar["Dict"])
    d2       = jsonDict["basic"]
    dGrid    = jsonDict["grid"]

    # Read all filters:
    f_eq = getFilter("filter_eq",jsonDict)
    f_gt = getFilter("filter_gt",jsonDict)
    f_lt = getFilter("filter_lt",jsonDict)

    # Get file naming pattern (needs YYYY MM and DD in there)
    fPattern = jsonDict["filePattern"]
    
    # Get main folder for files:
    folder   = jsonDict["folder"]

    NCDict = Dict{String, NCDatasets.CFVariable}()
    println("Input variables to output variables:")
    for (key, value) in dGrid
        println(key," : ", value)
        NCDict[key] = defVar(dsOut,key,Float32,("time","lat","lon"),deflatelevel=4, fillvalue=-9999)
        if ar["compSTD"]
            key2 = key*"_std"
            NCDict[key2] = defVar(dsOut,key2,Float32,("time","lat","lon"),deflatelevel=4, fillvalue=-9999, comment="Standard Deviation from data")
        end
    end
    println("")
    
    #dSIF = defVar(dsOut,"sif",Float32,("lat","lon"),deflatelevel=4, fillvalue=-9999)
    dN = defVar(dsOut,"n",Float32,("time","lat","lon"),deflatelevel=4, fillvalue=-9999, units="", long_name="Number of pixels in average")
    
    # Define data array
    mat_data          = zeros(Float32,(length(lat),length(lon),length(dGrid)))
    mat_data_variance = zeros(Float32,(length(lat),length(lon),length(dGrid)))
    mat_data_weights  = zeros(Float32,(length(lat),length(lon)))

    # Still hard-coded here, can be changed:
    nGrid = 10;

    points = zeros(Float32,(nGrid,nGrid,2))
    #global indices = zeros(Int32,(nGrid,nGrid,2))

    # Just to make sure we fill in attributes first time we read actual data:
    fillAttrib = true;

    # Loop through time:
    # Time counter
    cT = 1

    for chunk in 1:size(dateChunks, 1)
        files = String[];
        firstDate = dateChunks[chunk:chunk, :][1]
        lastDate  = dateChunks[chunk:chunk, :][2]
        println("Sub date range is: ", Date(firstDate), " to ", Date(lastDate))
        for day in firstDate:Dates.Day(1):lastDate
            filePattern = reduce(replace,["YYYY" => lpad(Dates.year(day),4,"0"), "MM" => lpad(Dates.month(day),2,"0"),  "DD" => lpad(Dates.day(day),2,"0")], init = fPattern)
            files       = [files;glob(filePattern, folder)]
        end
        
        fileSize = Int[];
        for f in files
            fileSize = [fileSize;stat(f).size]
        end

        # Loop through all files
        for a in files[fileSize.>0]

            fin = Dataset(a)
            # println("Read, ", a)
            
            # Read lat/lon bounds (required, maybe can change this to simple gridding in the future with just center):
            lat_in_ = getNC_var(fin, d2["lat_bnd"],true)
            lon_in_ = getNC_var(fin, d2["lon_bnd"],true)
            
            #println("Read")
            dim = size(lat_in_)

            # Transpose if orders are swapped
            if dim[1]==4
                lat_in_ = lat_in_'
                lon_in_ = lon_in_'
            end

            # Find all indices within lat/lon bounds:
            minLat = minimum(lat_in_, dims=2)
            maxLat = maximum(lat_in_, dims=2)
            minLon = minimum(lon_in_, dims=2)
            maxLon = maximum(lon_in_, dims=2)

            # Get indices within the lat/lon bounding box and check filter criteria (the last one filters out data crossing the date boundary):
            bool_add = (minLat[:,1].>latMin) .+ (maxLat[:,1].<latMax) .+ (minLon[:,1].>lonMin) .+ (maxLon[:,1].<lonMax) .+ ((maxLon[:,1].-minLon[:,1]).<50)
            
            bCounter = 5
            # Look for equalities
            for (key, value) in f_eq
                #println(key, " ", value)
                bool_add += (getNC_var(fin, key,false).==value)
                bCounter+=1
            end
            # Look for >
            for (key, value) in f_gt
                bool_add += (getNC_var(fin, key,false).>value)
                bCounter+=1
            end
            # Look for <
            for (key, value) in f_lt
                bool_add += (getNC_var(fin, key,false).<value)
                bCounter+=1
            end

            # If all were true, bool_add would be bCounter!
            idx = findall(bool_add.==bCounter)

            # Read data only for non-empty indices
            if length(idx) > 0
                #print(size(lat_in_))
                mat_in =  zeros(Float32,(length(lat_in_[:,1]),length(dGrid)))
                dim = size(mat_in)
                
                # Read in all entries defined in JSON file:
                co = 1
                
                # Do this only once:
                if fillAttrib
                    for (key, value) in dGrid
			# Need to change this soon to just go over keys(attrib), not this hard-coded thing. Just want to avoid another fill_value!
                        attribs = ["units","long_name","valid_range","description","unit","longname"]
                        for at in attribs
                            try
                                NCDict[key].attrib[at] = getNC_attrib(fin, value, at)
                            catch e
                                # @show e
                                # println(" Couldn't write attrib ", at)
                            end
                        end
                    end
                    fillAttrib=false
                end

                for (key, value) in dGrid
                    #println(key, value)
                    mat_in[:,co]=getNC_var(fin, value,false)
                    co += 1
                end

                iLat_ = ((lat_in_[idx,:].-latMin)/(latMax-latMin)*length(lat)).+1
                iLon_ = ((lon_in_[idx,:].-lonMin)/(lonMax-lonMin)*length(lon)).+1

                # @time  favg_all!(mat_data, mat_data_variance, mat_data_weights, ar["compSTD"], iLat_,iLon_,mat_in[idx,:],length(idx),dim[2],nGrid, latMin, latMax, lonMin,lonMax, length(lat), length(lon), points )
                favg_all!(mat_data, mat_data_variance, mat_data_weights, ar["compSTD"], iLat_,iLon_,mat_in[idx,:],length(idx),dim[2],nGrid, latMin, latMax, lonMin,lonMax, length(lat), length(lon), points )

                # println("Read ", a, " ", length(idx))
                println("Read ", basename(a))
            else
                # println("Read ", a, " ", length(idx))
                println("Read ", basename(a))
            end
            close(fin)
    #       catch
    #           println("Error in file caught")
    #       end
        end
        # Filter all data, set averages, still need to change row/column order here in the future!
        dims = size(mat_data)
        println("Averaging time step...")
        println("")
        if maximum(mat_data_weights) > 0
            dN[cT,:,:] = mat_data_weights
            dsTime[cT] = firstDate
            co = 1
            for (key, value) in dGrid
                da = round.(mat_data[:,:,co],sigdigits=6)
                da[mat_data_weights .< 1e-10] .= -9999
                NCDict[key][cT,:,:] = da
                if ar["compSTD"]
                    da = round.(sqrt.(mat_data_variance[:,:,co] ./ mat_data_weights), sigdigits = 6)
                    da[mat_data_weights.<1e-10] .= -9999
                    key2 = key * "_std"
                    NCDict[key2][cT,:,:] = da
                end
                #NCDict[key][cT,:,:]=da
                co += 1
            end
        else
            dN[cT,:,:] = 0
            dsTime[cT] = firstDate
        end
        cT += 1
        fill!(mat_data,0.0)
        fill!(mat_data_weights,0.0)
        fill!(mat_data_variance,0.0)
    end
    close(dsOut)

    # Permute?
    if ar["permute"] == true
        println("")
        println("--permute true. Permuting using nco in command line.")
        println("")
        nc_in = ar["outFile"]
        run(`ncpdq -a time,lat,lon -O $nc_in $nc_in`)
    end
    println("Gridding Finished! Time taken was: ")
    println(Dates.canonicalize(now() - startTime))
    println("")
end

main()
