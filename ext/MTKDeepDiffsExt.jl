module MTKDeepDiffsExt

using DeepDiffs, ModelingToolkit
using ModelingToolkit.BipartiteGraphs: Label,
    BipartiteAdjacencyList, unassigned,
    HighlightInt
using ModelingToolkit: SystemStructure,
    MatchedSystemStructure,
    SystemStructurePrintMatrix

"""
A utility struct for displaying the difference between two HighlightInts.

# Example
```julia
using ModelingToolkit, DeepDiffs

old_i = HighlightInt(1, :default, true)
new_i = HighlightInt(2, :default, false)
diff = HighlightIntDiff(new_i, old_i)

show(diff)
```
"""
struct HighlightIntDiff
    new::HighlightInt
    old::HighlightInt
end

function Base.show(io::IO, d::HighlightIntDiff)
    p_color = d.new.highlight
    (d.new.match && !d.old.match) && (p_color = :light_green)
    (!d.new.match && d.old.match) && (p_color = :light_red)

    (d.new.match || d.old.match) && printstyled(io, "(", color = p_color)
    if d.new.i != d.old.i
        Base.show(io, HighlightInt(d.old.i, :light_red, d.old.match))
        print(io, " ")
        Base.show(io, HighlightInt(d.new.i, :light_green, d.new.match))
    else
        Base.show(io, HighlightInt(d.new.i, d.new.highlight, false))
    end
    (d.new.match || d.old.match) && printstyled(io, ")", color = p_color)
end

"""
A utility struct for displaying the difference between two
BipartiteAdjacencyList's.

# Example
```julia
using ModelingToolkit, DeepDiffs

old = BipartiteAdjacencyList(...)
new = BipartiteAdjacencyList(...)
diff = BipartiteAdjacencyListDiff(new, old)

show(diff)
```
"""
struct BipartiteAdjacencyListDiff
    new::BipartiteAdjacencyList
    old::BipartiteAdjacencyList
end

function Base.show(io::IO, l::BipartiteAdjacencyListDiff)
    print(io,
        LabelDiff(Label(l.new.match === true ? "∫ " : ""),
            Label(l.old.match === true ? "∫ " : "")))
    (l.new.match !== true && l.old.match !== true) && print(io, "  ")

    new_nonempty = isnothing(l.new.u) ? nothing : !isempty(l.new.u)
    old_nonempty = isnothing(l.old.u) ? nothing : !isempty(l.old.u)
    if new_nonempty === true && old_nonempty === true
        if (!isempty(setdiff(l.new.highlight_u, l.new.u)) ||
            !isempty(setdiff(l.old.highlight_u, l.old.u)))
            throw(ArgumentError("The provided `highlight_u` must be a sub-graph of `u`."))
        end

        new_items = Dict(i => HighlightInt(i, :nothing, i === l.new.match) for i in l.new.u)
        old_items = Dict(i => HighlightInt(i, :nothing, i === l.old.match) for i in l.old.u)

        highlighted = union(map(intersect(l.new.u, l.old.u)) do i
                HighlightIntDiff(new_items[i], old_items[i])
            end,
            map(setdiff(l.new.u, l.old.u)) do i
                HighlightInt(new_items[i].i, :light_green,
                    new_items[i].match)
            end,
            map(setdiff(l.old.u, l.new.u)) do i
                HighlightInt(old_items[i].i, :light_red,
                    old_items[i].match)
            end)
        print(IOContext(io, :typeinfo => typeof(highlighted)), highlighted)
    elseif new_nonempty === true
        printstyled(io, map(l.new.u) do i
                HighlightInt(i, :nothing, i === l.new.match)
            end, color = :light_green)
    elseif old_nonempty === true
        printstyled(io, map(l.old.u) do i
                HighlightInt(i, :nothing, i === l.old.match)
            end, color = :light_red)
    elseif old_nonempty !== nothing || new_nonempty !== nothing
        print(io,
            LabelDiff(Label(new_nonempty === false ? "∅" : "", :light_black),
                Label(old_nonempty === false ? "∅" : "", :light_black)))
    else
        printstyled(io, '⋅', color = :light_black)
    end
end

"""
A utility struct for displaying the difference between two Labels
in git-style red/green highlighting.

# Example
```julia
using ModelingToolkit, DeepDiffs

old = Label("before")
new = Label("after")
diff = LabelDiff(new, old)

show(diff)
```
"""
struct LabelDiff
    new::Label
    old::Label
end
function Base.show(io::IO, l::LabelDiff)
    if l.new != l.old
        printstyled(io, l.old.s, color = :light_red)
        length(l.new.s) != 0 && length(l.old.s) != 0 && print(io, " ")
        printstyled(io, l.new.s, color = :light_green)
    else
        print(io, l.new)
    end
end

"""
A utility struct for displaying the difference between two
(Matched)SystemStructure's in git-style red/green highlighting.

# Example
```julia
using ModelingToolkit, DeepDiffs

old = SystemStructurePrintMatrix(...)
new = SystemStructurePrintMatrix(...)
diff = SystemStructureDiffPrintMatrix(new, old)

show(diff)
```
"""
struct SystemStructureDiffPrintMatrix <:
       AbstractMatrix{Union{LabelDiff, BipartiteAdjacencyListDiff}}
    new::SystemStructurePrintMatrix
    old::SystemStructurePrintMatrix
end

function Base.size(ssdpm::SystemStructureDiffPrintMatrix)
    max.(Base.size(ssdpm.new), Base.size(ssdpm.old))
end

function Base.getindex(ssdpm::SystemStructureDiffPrintMatrix, i::Integer, j::Integer)
    checkbounds(ssdpm, i, j)
    if i > 1 && (j == 4 || j == 9)
        old = new = BipartiteAdjacencyList(nothing, nothing, unassigned)
        (i <= size(ssdpm.new, 1)) && (new = ssdpm.new[i, j])
        (i <= size(ssdpm.old, 1)) && (old = ssdpm.old[i, j])
        BipartiteAdjacencyListDiff(new, old)
    else
        old = new = Label("")
        (i <= size(ssdpm.new, 1)) && (new = ssdpm.new[i, j])
        (i <= size(ssdpm.old, 1)) && (old = ssdpm.old[i, j])
        LabelDiff(new, old)
    end
end

function DeepDiffs.deepdiff(old::Union{MatchedSystemStructure, SystemStructure},
        new::Union{MatchedSystemStructure, SystemStructure})
    new_sspm = SystemStructurePrintMatrix(new)
    old_sspm = SystemStructurePrintMatrix(old)
    Base.print_matrix(stdout, SystemStructureDiffPrintMatrix(new_sspm, old_sspm))
end

end # module
