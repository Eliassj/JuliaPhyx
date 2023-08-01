module JuliaPhyx

using DataFrames
const dt = DataFrames
using Gtk
using MAT
using NPZ
using CairoMakie
const plt = CairoMakie
using CSV
using Mmap
const mp = Mmap


# Initialize a PhyOutput object
struct PhyOutput
    _spiketimes::Matrix{UInt64}
    _info::DataFrames.DataFrame
    _meta::Dict{SubString{String}, SubString{String}}
    _binpath::String

    """
    PhyOutput(
        phydir::String = "",
        glxdir::String = "",
        same::Bool = false
    )


`phydir`: Optional. The directory containing kilosort/phy output.
`glxdir`: Optional. The directory containing spikeGLX output (.ap.bin/.ap.meta)
`same`: If true spikeGLX output will be assumed to be in the same directory as Kilosort/Phy output.

Create a PhyOutput struct containing spiketimes(_spiketimes), info(_info), spikeGLX metadata(_meta) and the path to the associated .ap.bin file(_binpath).
"""
    function PhyOutput(
        phydir::String = "",
        glxdir::String = "",
        same::Bool = false
    )
        if length(phydir) == 0
            println("Select phy/kilosort output directory")
            phydir = Gtk.open_dialog_native("Select Kilosort/Phy output Folder", action=GtkFileChooserAction.SELECT_FOLDER)
        end
        if length(glxdir) == 0 && same == false
            println("Select spikeGLX output directory")
            glxdir = Gtk.open_dialog_native("Select spikeGLX output directory", action=GtkFileChooserAction.SELECT_FOLDER)
        end
        if length(glxdir) == 0 && same == true
            glxdir = phydir
        end
        println("Importing good clusters")
        clusters::Vector{UInt64} = convert(Vector{UInt64}, NPZ.npzread(phydir*"\\spike_clusters.npy"))
        times::Vector{UInt64} = NPZ.npzread(phydir*"\\spike_times.npy")[:,1]
        spiketimes::Matrix{UInt64} = [clusters times]
        info::DataFrames.DataFrame = CSV.read(phydir*"\\cluster_info.tsv", DataFrame)

        isgood(group) = group == "good"
        info = subset(info, :group => ByRow(isgood), skipmissing = true)

        ininfo(cluster) = cluster in info[!, 1]
        spiketimes = spiketimes[ininfo.(spiketimes[:,1]),:]

        glxfiles = readdir(glxdir, join = true)
        binfile::String = [f for f in glxfiles if f[length(f)-6:length(f)] == ".ap.bin"][1]
        metafile::String = [f for f in glxfiles if f[length(f)-7:length(f)] == ".ap.meta"][1]

        # Read metadata
        tmp = open(metafile, "r")
        metaraw = readlines(tmp)
        close(tmp)
        metaraw = split.(metaraw, "=")
        metadict = Dict(i[1] => i[2] for i in metaraw)

        new(spiketimes, info, metadict, binfile)

        end
end




function spikemmap(p::PhyOutput)
    n::Int = parse(Int, p._meta["nSavedChans"])
    s::Int = Int(parse(Int, p._meta["fileSizeBytes"]) / (2*n))
    tmp::IOStream = open(p._binpath, "r")
    m::Array{Int16, 2} = mp.mmap(tmp, Array{Int16, 2}, (n, s), 0)
    close(tmp)
    return m
end

function tovolts(meta::Dict{SubString{String}, SubString{String}}, i::Union{Vector{Int16}, Array{Int16}})
    Imax::Float64 = parse(Float64, meta["imMaxInt"])
    Vmax::Float64 = parse(Float64, meta["imAiRangeMax"])
    if meta["imDatPrb_type"] == "0"
        tbl = meta["~imroTbl"]
        tbl = split(tbl, ")(")
        tbl = split(tbl[2], " ")
        gain = parse(Int, tbl[4])
    else
        gain::Float64 = 80
    end
    cfactor::Float64 = Vmax / Imax / gain
    return i .*cfactor
end

"""
    getchan(
    p::PhyOutput,
    ch::Union{Int, Vector{Int}, UnitRange{Int64}},
    tmin::Union{Float64, Int},
    tmax::Union{Float64, Int},
    converttoV::Bool = true
    )

`p`: Created by JuliaPhyx.PhyOutput()
`ch`: Channel index (1-based) to extract. May be a single `Int`, a `Vector{Int}` or a `UnitRange{Int64}`.
`tmin`: Time in seconds to begin from.
`tmax`: Time in seconds to stop at.
`ConverttoV`: Should the values be converted to volts?

Uses the path stored in p._binfile and p._meta to import channels from the correct binfile in the correct way.

See also [`PhyOutput`](@ref)
"""
function getchan(
    p::PhyOutput,
    ch::Union{Int, Vector{Int}, UnitRange{Int64}},
    tmin::Union{Float64, Int},
    tmax::Union{Float64, Int, String},
    converttoV::Bool = true
    )

    filemax::Int64 = Int64(parse(Int64, p._meta["fileSizeBytes"]) / parse(Int64, p._meta["nSavedChans"]) / 2)
    samplfrq::Float64 = parse(Float64, p._meta["imSampRate"])
    # Convert times (IN SECONDS) to sample frequencies
    # Add 1 since indexing is 1-based(?)

    #Check if tmin < 0
    if tmin < 0
        ArgumentError("tmin must be above 0")
    end

    # Check if max time should be used
    if typeof(tmax) == String
        if tmax != "max"
            ArgumentError("tmax must be a number or the string 'max'")
        end
        tmax = filemax
    else
        # Convert tmax to sample frequency
        tmax::Int64 = Int(round((tmax * samplfrq))) + 1::Int64
    end
    # Check if tmin >= tmax
    if tmin >= tmax
        ArgumentError("tmin must be <tmax")
    end

    # Convert tmin to sample frequency
    tmin::Int64 = Int(round((tmin * samplfrq))) + 1::Int64
    
    if tmax > filemax
        ArgumentError("tmax larger than length of recording")
    end

    # Memory map data and load the selected chunks
    karta::Matrix{Int16} = spikemmap(p)
    len::UnitRange{Int64} = tmin:tmax

    r::Union{Matrix{Int16}, Vector{Int16}} = karta[ch, len]

    if converttoV
        conv::Union{Matrix{Float64}, Vector{Float64}} = tovolts(p._meta, r)
        return conv
    else
        return r
    end

    end

    function writechanbin(
        channel::Union{Matrix{Float64}, Matrix{Int16}, Vector{Int16}, Vector{Float64}},
        dest::String="")

        if dest == ""
            dest = Gtk.open_dialog_native("Select save location", action=GtkFileChooserAction.SAVE)
        end
        if dest[end-3:end] != ".bin"
            dest = dest*".bin"
        end
        open(dest, "w") do file
            write(file, channel)            
            close(file)
        end   
        println("written to: ", dest)
    end

    function importchanint16(path::String="", format = "Int16" )
        if path == ""
            path = Gtk.open_dialog_native("Select trigger file", action=GtkFileChooserAction.GTK_FILE_CHOOSER_ACTION_OPEN)
        end
        if path[end-3:end] == ".bin"
            tmp::IOStream = open(path, "r")
            res = reinterpret(Int16, read(tmp))
            close(tmp)
            return res
        elseif path[end-3:end] == ".csv"
            res2::Matrix = Matrix(CSV.read(path, DataFrame, header = false))
            return res2
        else
            ArgumentError("Only '.bin' or '.csv' files allowed.")
        end
    end

end