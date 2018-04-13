

function filefromid(ID::Integer,DataDir)
    tarr = []
    filearr = readdir("$DataDir") 
    for i in filearr
        fi = stat("$DataDir/$i")
        datei = Dates.unix2datetime(fi.mtime)
        hh = parse(Int,Dates.format(datei,"HH"))+1
        ti = parse(Int,Dates.format(datei,"yymmdd")*"$hh"*Dates.format(datei,"MM"))
        #println("$i : $ti")
        push!(tarr,ti)
    end

    i = findfirst(tarr,ID)

    if i==0
        error("No sample found")
    end

    FNAME = "$DataDir/$(filearr[i])"

    return FNAME
end

# ID = 1801161626
# FNAME = filefromid(ID)
