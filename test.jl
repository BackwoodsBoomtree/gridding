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

folder   = "/mnt/g/TROPOMI/caltech/original"
fPattern = "TROPO_SIF_YYYY-MM-DD_*.nc"

startDate = DateTime("2019-12-19")
stopDate  = DateTime("2021-01-11")
dDay      = Dates.Day(8)

modDateList = collect(Date(startDate):Dates.Day(8):Date(stopDate))


startDate = "2019-12-19"
stopDate  = "2020-01-11"

function modDates(startDate, endDate)
    startDate = DateTime(startDate)
    endDate   = DateTime(endDate)
    nyear     = Dates.Year(endDate).value - Dates.Year(startDate).value + 1
    yearList  = [Dates.Year(startDate).value : 1 : Dates.Year(endDate).value;]
    println(yearList)
    println("Number of years is ", nyear)
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
        dateList = replace(dateList, DateTime(string(year, "-01-02")) => DateTime(string(year, "-12-31")))
        dateList = replace(dateList, DateTime(string(year, "-01-03")) => DateTime(string(year, "-12-31")))
    end
    return(dateList)
end

myList = modDates("2019-12-19", "2021-01-11")


first    = DateTime("2019-12-19");
last     = DateTime("2021-01-31");

firstList = collect(first:Dates.Day(8):last)
endList   = firstList + Dates.Day(7)
dateList = hcat(firstList, endList)

# If modLike but start date is bad then exit,
# else calculate cT (timesteps) length
if ar["modLike"]
    modstartYear = Dates.Year(startDate).value
    modDateList  = collect(Date(string(modstartYear, "-01-01")):Dates.Day(8):Date(string(modstartYear, "-12-31")))
    if Date(startDate) ∉ modDateList
        println("ERROR: Input --startDate does not fall on an 8-day MODIS start date.")
        exit()
    else
        nyear     = Dates.Year(stopDate).value - Dates.Year(startDate).value + 1
        if nyear > 1            
            # Length of first year
            cT = length(collect(Date(startDate):Dates.Day(8):Date(string(Dates.Year(startDate).value, "-12-31"))))
            # Add Length of last year
            cT = cT + length(collect(Date(string(Dates.Year(stopDate).value, "-01-01")):Dates.Day(8):Date(stopDate)))
            # Add length of remaining years (always 46)
            if nyear > 2
                cT = cT + 46 * (nyear - 2)
            end
        end
        println("Temporal resolution is modis-like. Number of time steps is ", cT)
    end
end




startYear   = Dates.Year(startDate).value
modDateList = collect(Date(string(2021, "-01-01")):Dates.Day(8):Date(string(2021, "-04-01")))

modDateList = collect(Date(startDate):Dates.Day(8):Date(stopDate))

nyear    = Dates.Year(stopDate).value - Dates.Year(startDate).value + 1









if nyear > 1
    # Length of first year
    cT = length(collect(Date(startDate):Dates.Day(8):Date(string(Dates.Year(startDate).value, "-12-31"))))
    # Length of last year
    cT = cT + length(collect(Date(string(Dates.Year(stopDate).value, "-01-01")):Dates.Day(8):Date(stopDate)))
    if nyear > 2
        cT = cT + 46 * (nyear - 2)
    end
end

if 1 == 2 && 2 == 2
    println("okay")
else
    println("not okay")
end



if Dates.Year(startDate).value < Dates.Year(stopDate).value
    # Create list of years
    nyear    = Dates.Year(stopDate).value - Dates.Year(startDate).value + 1
    yearList = Int64[];
    for i in 1:nyear
        push!(yearList, Dates.Year(startDate).value + i - 1)
    end

    for y in 1:length(yearList)
        if y == 1
            yearStartDate = startDate
            yearStopDate  = DateTime(string(yearList[y], "-12-31"))

    end
end


1+1

a = 1
b = 2
c = 11

if 1 == 1
    a = 2
    print(a)
end

print(a)

for i in a:b:c 
    println(i) 
    println(a)
    # a = a + 1
end




for d in startDate:dDay:stopDate
    files = String[];
    for di in d:Dates.Day(1):d + dDay - Dates.Day(1)
        if di <= stopDate
            filePattern = reduce(replace,["YYYY" => lpad(Dates.year(di),4,"0"), "MM" => lpad(Dates.month(di),2,"0"),  "DD" => lpad(Dates.day(di),2,"0")], init=fPattern)
            println(filePattern)
            files = [files;glob(filePattern, folder)]
        end
    end
    fileSize = Int[];
    for f in files
        fileSize = [fileSize;stat(f).size]
    end
    if length(files) < Dates.value(dDay)
        println("Warning: ", length(files), " files used for gidding for date range ", Date(d), " to ", Date(d + dDay - Dates.Day(1)))
        println("First and last files used were: ", basename(files[1]), " to ", basename(last(files)))
        println("")
    end
end

# cT = 1
# for d in startDate:dDay:stopDate
#     files = String[];
#     for di in d:Dates.Day(1):d+dDay-Dates.Day(1)
#         #println("$(@sprintf("%04i-%02i-%02i", Dates.year(di),Dates.month(di),Dates.day(di)))")

#         filePattern = reduce(replace,["YYYY" => lpad(Dates.year(di),4,"0"), "MM" => lpad(Dates.month(di),2,"0"),  "DD" => lpad(Dates.day(di),2,"0")], init=fPattern)
#         #println(filePattern, " ", folder)
#         files = [files;glob(filePattern, folder)]
#     end
#     fileSize = Int[];
#     for f in files
#         fileSize = [fileSize;stat(f).size]
#     end
#     println(files)