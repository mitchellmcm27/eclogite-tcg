"""
MIT License

Copyright (c) 2018-2021 C. Brenhin Keller and contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

using StatGeochem
using DelimitedFiles

const Collection{T} = Union{AbstractArray{<:T}, NTuple{N,T}} where N
"""

"""
Set up a PerpleX calculation for a single bulk composition across an entire
2d P-T space. P specified in bar and T in Kelvin.
"""
function perplex_build_vertex(perplexdir::String, scratchdir::String, composition::Collection{Number},
        elements::Collection{String},
        P::NTuple{2,Number}=(280, 28000), T::NTuple{2,Number}=(273.15, 1500+273.15);
        dataset::String="hp11ver.dat",
        xnodes::Integer=40,
        ynodes::Integer=40,
        solution_phases::String="O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\ncAmph(DP)\nT\nB\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar_B\nDo(HP)\nF\n",
        excludes::String="ts\nparg\ngl\nged\nfanth\ng\n",
        mode_basis::String="vol",
        composition_basis::String="wt",
        fluid_eos::Number=5,
        saturated_fluid::String="",
        name::String="scratch"
    )

    xnodes_r = xnodes .* 4
    ynodes_r = ynodes .* 4
    bymass = "y"
    if(composition_basis=="mol")
        bymass = "n"
    end
    # Use chemical potentials, activities or fugacities as independent variables (Y/N)?
    chempot = "n"

    build = joinpath(perplexdir, "build")# path to PerpleX build
    vertex = joinpath(perplexdir, "vertex")# path to PerpleX vertex

    #Configure working directory
    prefix = joinpath(scratchdir, "$(name)/")
    system("rm -rf $prefix; mkdir -p $prefix")

    # Place required data files
    system("cp $(joinpath(perplexdir,dataset)) $prefix")
    system("cp $(joinpath(perplexdir,"perplex_option.dat")) $prefix")
    if (dataset=="stx21ver.dat")
        system("cp $(joinpath(perplexdir,"stx21_solution_model.dat")) $prefix/solution_model.dat")
    else
        system("cp $(joinpath(perplexdir,"solution_model.dat")) $prefix")
    end

    # Edit data files to specify number of nodes at which to solve
    system("sed -e \"s/x_nodes .*|/x_nodes                   $xnodes $xnodes_r |/\" -i.backup $(prefix)perplex_option.dat")
    system("sed -e \"s/y_nodes .*|/y_nodes                   $ynodes $ynodes_r |/\" -i.backup $(prefix)perplex_option.dat")

    # Specify whether we want volume or weight percentages
    system("sed -e \"s/proportions .*|/proportions                    $mode_basis |/\" -i.backup $(prefix)perplex_option.dat")
    #system("sed -e \"ss/composition_system .*|/composition_system             $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")
    #system("sed -e \"s/composition_phase .*|/composition_phase              $composition_basis |/\" -i.backup $(prefix)perplex_option.dat")

    # Create build batch file
    # Options based on Perplex v6.8.7
    fp = open(prefix*"build.bat", "w")

    # Name, components, and basic options. P-T conditions.
    # default fluid_eos = 5: Holland and Powell (1998) "CORK" fluid equation of state
    elementstring = join(elements .* "\n")
    write(fp,"$name\n$dataset\nperplex_option.dat\nn\n2\n$saturated_fluid\nn\nn\n$elementstring\nn\n2\n$(first(T))\n$(last(T))\n$(first(P))\n$(last(P))\n$(bymass)\n") # v6.8.7

    # Whole-rock composition
    for i ∈ eachindex(composition)
        write(fp,"$(composition[i]) ")
    end
    # Solution model
    printfile = "n"
    excl_endmember = "y"
    prompt_phases = "n"
    incl_soln = "y"
    write(fp,"\n$printfile\n$excl_endmember\n$prompt_phases\n$excludes\n$incl_soln\nsolution_model.dat\n$solution_phases\n$(name)_Pseudosection")
    close(fp)

    # build PerpleX problem definition
    system("cd $prefix; $build < build.bat > build.log")

    # Run PerpleX vertex calculations
    result = system("cd $prefix; echo $name | $vertex > vertex.log")
    return result
end

export perplex_build_vertex

function perplex_pssect(perplexdir::String, scratchdir::String;
    name::String="scratch")
    # Run PSSECT to generate a pseudosection

    pssect = joinpath(perplexdir, "pssect")# path to PerpleX werami
    prefix = joinpath(scratchdir, "$(name)/") # path to data files

    # Create werami batch file
    fp = open(prefix*"pssect.bat", "w")
    # v6.7.8 pseudosection
    write(fp,"$name\nn\n")
    close(fp)

    # Make sure there isn"t already an output
    system("rm -f $(prefix)$(name).ps")

    # Run PSSECT
    system("cd $prefix; $pssect < pssect.bat > pssect.log")
end

export perplex_pssect

function perplex_werami_profile(perplexdir::String, scratchdir::String, P::NTuple{2,Number}, T::NTuple{2,Number};
    name::String="scratch", npoints::Integer=100, include_fluid="y", importas=:Dict)
    # Query a new path from a pseudosection

    werami = joinpath(perplexdir, "werami")# path to PerpleX werami
    prefix = joinpath(scratchdir, "$(name)/") # path to data files

    # Create werami batch file
    fp = open(prefix*"werami.bat", "w")
    # v6.7.8 pseudosection
    write(fp,"$name\n3\nn\n$(first(T))\n$(first(P))\n$(last(T))\n$(last(P))\n$npoints\n25\nn\n$include_fluid\n0\n")
    close(fp)

    # Make sure there isn"t already an output
    system("rm -f $(prefix)$(name)_1.tab*")

    # Extract Perplex results with werami
    system("cd $prefix; $werami < werami.bat > werami.log")

    # Ignore initial and trailing whitespace
    system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(name)_1.tab")
    # Merge delimiters
    system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(name)_1.tab")

    # Read results and return them if possible
    result = importas==:Dict ? Dict() : ()
    try
        # Read data as an Array{Any}
        data = readdlm("$(prefix)$(name)_1.tab", ' ', skipstart=8)
        # Convert to a dictionary
        result = elementify(data, sumduplicates=false, importas=importas)
    catch e
        showerror(stdout, e)
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(name)_1.tab could not be parsed, perplex may not have run"
    end
    return result
end

export perplex_werami_profile

function perplex_werami_rho(perplexdir::String, scratchdir::String;
    include_fluid::String="", name::String="scratch", importas=:Dict)
    # Query a new path from a pseudosection

    werami = joinpath(perplexdir, "werami")# path to PerpleX werami
    prefix = joinpath(scratchdir, "$(name)/") # path to data files

    # Create werami batch file
    fp = open(prefix*"werami.bat", "w")
    # v6.7.8 pseudosection
    write(fp,"$name\n2\n2\nn\n$include_fluid\n0\nn\n1\n0\n")
    close(fp)

    # Make sure there isn"t already an output
    system("rm -f $(prefix)$(name)_2.tab*")

    # Extract Perplex results with werami
    system("cd $prefix; $werami < werami.bat > werami.log")

    # Ignore initial and trailing whitespace
    system("sed -e \"s/^  *//\" -e \"s/  *\$//\" -i.backup $(prefix)$(name)_2.tab")
    # Merge delimiters
    system("sed -e \"s/  */ /g\" -i.backup $(prefix)$(name)_2.tab")

    # Read results and return them if possible
    result = importas==:Dict ? Dict() : ()
    try
        # Read data as an Array{Any}
        data = readdlm("$(prefix)$(name)_2.tab", ' ', skipstart=12)
        # Convert to a dictionary
        result = elementify(data, sumduplicates=false, importas=importas)
    catch e
        showerror(stdout, e)
        # Return empty dictionary if file doesn't exist
        @warn "$(prefix)$(name)_2.tab could not be parsed, perplex may not have run"
    end
    return result
end

export perplex_werami_rho

function perplex_werami_point(perplexdir::String, scratchdir::String, P::Number, T::Number; name::String="scratch")
    werami = joinpath(perplexdir, "werami")# path to PerpleX werami
    prefix = joinpath(scratchdir, "$(name)/") # path to data files

    # Sanitize T inputs to avoid PerpleX escape sequence
    if P == 99
        P = 99.001
    end
    if T == 99
        T = 99.001
    end

    # Create werami batch file
    # Options based on Perplex v6.7.2
    fp = open(prefix*"werami.bat", "w")
    write(fp,"$name\n1\n$T\n$P\n99\n99\n0\n")
    close(fp)

    # Make sure there isn"t already an output
    system("rm -f $(prefix)$(name)_1.txt")

    # Extract Perplex results with werami
    system("cd $prefix; $werami < werami.bat > werami.log")

    # Read results and return them if possible
    data = ""
    try
        # Read entire output file as a string
        fp = open("$(prefix)$(name)_1.txt", "r")
        data = read(fp, String)
        close(fp)
    catch e
        showerror(stdout, e)
        # Return empty string if file doesn't exist
        @warn "$(prefix)$(name)_1.txt could not be parsed, perplex may not have run"

    end
    return data
end
export perplex_werami_point