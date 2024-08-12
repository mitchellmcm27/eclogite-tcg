using StatGeochem
using JSON

include("perplex_api.jl")

# Absolute path to Perple_X installation
resourcepath = "/root/resources"
perplexdir = joinpath(resourcepath,"perplex-stable")

mode_basis = "vol"

T_range_2d = (300+273.15, 1300+273.15) # Kelvin
P_range_2d = (5000, 30000) # bar

T_range_1d = (650+273.15, 850+273.15) # Kelvin
P_range_1d = (5000, 30000) # bar

T_point = 900+273.15 # K
P_point = 2.0e4 # bar

T_surf = 273.15

force_pseudosection = false

args = ARGS
if "-f" in ARGS
    force_pseudosection = true
    args = filter(x->x!="-f", args)
end

compositions = JSON.parsefile("compositions.json")
all_comp_names = collect(keys(compositions))
comp_names = []
if length(args)>0
    for arg in args
        global comp_names = vcat(comp_names, arg)
    end
else 
    global comp_names = all_comp_names
end

for name in comp_names
    composition_name = name
    comp = compositions[name]
    dataset = comp["dataset"] isa String ? comp["dataset"] : "stx21ver"
    scratchdir = "./output/"
    if dataset == "stx21ver"
        if occursin("pyrolite",name) || occursin("harzburgite",name)
            println("excluding no phases")
            phases_exclude = []
        else
            println("excluding majoritic garnet endmembers")
            phases_exclude = ["maj","namaj","namj"]
        end
        # all solution phases
        phases = ["O","Pl","Sp","Cpx","Wad","Ring","Pv","Wus","C2/c","Opx","Aki","Ppv","CF","Gt","NaAl"]
        include_fluid="" # perplex will not ask, so leave it blank
        saturated_fluid=""
    else
        # Follows Li et al. 2022
        # Not used for Stixrude database
        phases_exclude = [
            #melt(G)
            "sil8L",
            "ctjL",
            "fo8L",
            "wi8L",
            "fa8L",
            "q8L",
            "qjL",
            "dijL",
            "jdjL",
            "ctjL",
            "fojL",
            "fajL",
            "hmjL",
            "ekjL",
            #cAmph(G)
            #Omph(HP)
            #Gt(W)
            "andr",
            #Opx(W)
            #Bi(W)
            #Mica(W)
            "ma",
            #Ilm(WPH)
            "ilm",
            "ilm_nol",
            #T
            #Sp(WPC)
            #Ep(HP11)
            #Pl(I1,HP)
            #Fsp(C1)
            ]
        phases = [
            "melt(G)",
            "cAmph(G)",
            "Omph(HP)",
            "Gt(W)",
            "Opx(W)" ,
            "Bi(W)",
            "Mica(W)" ,
            "Ilm(WPH)", # ilmenite-hematite
            "T",
            "Sp(WPC)", # magnetite-spinel
            "Ep(HP11)",
            "Pl(I1,HP)",
            "Fsp(C1)"
        ]
        include_fluid="n" # perplex asks
        saturated_fluid="n"
    end
    oxide_comp = convert(Array{Number},comp["composition"])
    oxides = convert(Array{String},comp["elements"])
    if dataset == "stx21ver"
        oxides = map((s) -> uppercase(s), oxides)
    end
    composition_basis = comp["basis"]

    println(sum(oxide_comp))
    println(oxide_comp)
    println(phases)
    println(oxides)
    println(oxide_comp)

    blkpath = joinpath(scratchdir,composition_name,composition_name*".blk")
    pseudosection_exists = isfile(blkpath)

    excludes = join(phases_exclude,"\n")
    if (length(phases_exclude)>0)
        excludes *= "\n"
    end
    solution_phases = join(phases,"\n")*"\n"

    if force_pseudosection || !pseudosection_exists
        print("Solving pseudosection...\n")
        perplex_build_vertex(perplexdir, scratchdir, oxide_comp,oxides, P_range_2d, T_range_2d, 
            dataset=dataset*".dat", 
            excludes=excludes,
            solution_phases=solution_phases, 
            composition_basis=composition_basis, 
            mode_basis=mode_basis, 
            name=composition_name,
            saturated_fluid=saturated_fluid
        )
        perplex_pssect(perplexdir, scratchdir, name=composition_name)
    end

    print("Extracting profile...\n")
    perplex_werami_profile(perplexdir, scratchdir, P_range_1d, T_range_1d,name=composition_name)

    print("Extracting density grid...\n")
    perplex_werami_rho(perplexdir, scratchdir, include_fluid=include_fluid, name=composition_name)

    print("Extracting point...\n")
    perplex_werami_point(perplexdir,scratchdir,P_point,T_point,name=composition_name)
    print(point)
end