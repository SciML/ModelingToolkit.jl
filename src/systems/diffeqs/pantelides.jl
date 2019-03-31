struct BiGraph{TA,TB,TC}
    eqs::TA
    vars::TB
    edges::TC
end

function init_assign(G)
    assign = Dict{eltype(G.vars), Union{Nothing, eltype(G.eqs)}}()
    # second function
    for j in G.vars
        assign[j] = nothing
    end
    return assign
end

function augmentpath!(i, pathfound, color, assign, vars, edges)
    union!(color, i)
    idx = -1
    for j in vars
        assign[j] === nothing || continue
        (i => j) in edges && (idx = j; break)
    end
    if idx != -1
        pathfound[] = true
        assign[idx] = i
        return nothing
    end
    for j in vars
        !((i => j) in edges && !(j in color)) && continue
        union!(color, j)
        k = assign[j]
        augmentpath!(k, pathfound, color, assign, vars, edges)
        pathfound[] && (assign[j] = i; return nothing)
    end
    return nothing
end

function construct_augmentpath!(G, assign, pathfound)
    #main operation in function
    color = Set(similar(G.eqs, 0))
    for i in G.eqs
        empty!(color)
        pathfound[] = false
        while !pathfound[]
            pathfound[] = false
            augmentpath!(i, pathfound, color, assign, G.vars, G.edges)
            pathfound[] || break
        end
    end
    return nothing
end

#function pantelides(m,n,G,diff)
#    assign = Dict()
#    b = Dict()
#    for j=1:m  # dict and vector() which is better. for loop?
#        assign[j] = 0
#    end
#    for j=1:n
#        b[j] = 0
#    end
#    pathfound = false
#    for k = 1:n
#        while !pathfound
#            # 3b-1
#            filter!(a -> a!=0,diff)
#            filter!(a -> a[2] in diff,G.edges)
#            color=Vector()
#            pathfound=Ref(false)
#            augmentpath!(k,pathfound,color,assign,G.vars,G.edges)
#            # 3b-5 to do
#            pathfound[] || return assign
#        end
#    end
#    return assign
#end
