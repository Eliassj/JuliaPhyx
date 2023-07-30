module JuliaPhyx

using DataFrames
const dt = DataFrames
using Gtk
using MAT
using NPZ
using Plots
const plt = Plots
using CSV
using Mmap
const mp = Mmap


# Initialize a PhyOutput object
struct PhyOutput
    _spiketimes
    _info
    _meta
    _binpath

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
        clusters = convert(Vector{UInt64}, NPZ.npzread(phydir*"\\spike_clusters.npy"))
        times = NPZ.npzread(phydir*"\\spike_times.npy")[:,1]
        spiketimes = [clusters times]
        info = CSV.read(phydir*"\\cluster_info.tsv", DataFrame)

        isgood(group) = group == "good"
        info = subset(info, :group => ByRow(isgood), skipmissing = true)

        ininfo(cluster) = cluster in info[!, 1]
        spiketimes = spiketimes[ininfo.(spiketimes[:,1]),:]

        glxfiles = readdir(glxdir, join = true)
        binfile = [f for f in glxfiles if f[length(f)-6:length(f)] == ".ap.bin"][1]
        metafile = [f for f in glxfiles if f[length(f)-7:length(f)] == ".ap.meta"][1]

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
    n = parse(Int, p._meta["nSavedChans"])
    s = Int(parse(Int, p._meta["fileSizeBytes"]) / (2*n))
    tmp = open(p._binpath, "r")
    m::Array{Int16, 2} = mp.mmap(tmp, Array{Int16, 2}, (n, s), 0)
    close(tmp)
    return m
end

function tovolts(meta::Dict{SubString{String}, SubString{String}}, i::Union{Vector{Int16}, Array{Int16}})
    Imax::Float64 = parse(Float64, meta["imMaxInt"])
    Vmax::Float64 = parse(Float64, meta["imAiRangeMax"])
    if meta["imDatPrb_type"] == "0"
        tbl=meta["~imroTbl"]
        tbl = split(tbl, ")(")
        tbl = split(tbl[2], " ")
        gain = parse(Int, tbl[4])
    else
        gain = 80
    end
    cfactor = Vmax / Imax / gain
    return i .*cfactor
end

function getchan(
    p::PhyOutput,
    ch::Union{Int, Vector{Int}, UnitRange{Int64}},
    tmin::Union{Float64, Int},
    tmax::Union{Float64, Int},
    converttoV::Bool = true
    )

    # Convert times (IN SECONDS) to sample frequencies
    # Add 1 since indexing is 1-based
    samplfrq::Float64 = parse(Float64, p._meta["imSampRate"])
    tmin::Int64 = Int(round((tmin * samplfrq) + 1))
    tmax::Int64 = Int(round((tmax * samplfrq) + 1))
    filemax::Int64 = Int64(parse(Int64, p._meta["fileSizeBytes"]) / parse(Int64, p._meta["nSavedChans"]) / 2)
    if tmax > filemax
        ArgumentError("tmax may at most be $(p._meta["fileTimeSecs"])")
    end
    # Memory map data and load the selected chunks
    karta = spikemmap(p)
    
    r = karta[ch, tmin:tmax]

    r = tovolts(p._meta, r)
    return r
end

end

Plots.pl