struct BiGraph
    nodese
    nodesv
    edges
end

function augment_path(G)
    assign = Dict()
    # second function
    function check_alg!(G,assign,pathfound)
        function augmentpath!(i,pathfound,color,assign,nodesv,edges)
            function checkcond(assign,nodesv,edges)
                for j in nodesv
                    if assign[j] == 0 && (count(x -> x[2]==j && x[1]==i,edges) > 0)
                        return j
                    end
                end
                return false
            end
            union!(color,i)
            j = checkcond(assign,nodesv,edges)
            if j != false
                pathfound[1]=true
                assign[j] = i
                return nothing
            end
            for j in nodesv
                if count(x -> x[2]==j && x[1]==i,edges) > 0 && !(j in color)
                    union!(color,j)
                    k=assign[j]
                    augmentpath!(k,pathfound,color,assign,nodesv,edges)
                    if pathfound[1]
                        assign[j] = i
                        return nothing
                    end
                end
            end
            return assign
        end
        #main operation in function
        for i in G.nodese
            pathfound[1] = false
            while !pathfound[1]
                color=Vector()
                pathfound[1] = false
                augmentpath!(i,pathfound,color,assign,G.nodesv,G.edges)
                pathfound[1] || break
            end
        end
    end
    for j in G.nodesv
        assign[j] = 0
    end
    check_alg!(G,assign,[false])
    return assign
end

function pantelides(m,n,G,diff)
    assign = Dict()
    b = Dict()
    for j=1:m  # dict and vector() which is better. for loop?
        assign[j] = 0
    end
    for j=1:n
        b[j] = 0
    end
    pathfound = false
    for k = 1:n
        while !pathfound
            # 3b-1
            filter!(a -> a!=0,diff)
            filter!(a -> a[2] in diff,G.edges)
            color=Vector()
            pathfound=false
            augmentpath!(k,pathfound,color,assign,G.nodesv,G.edges)
            # 3b-5 to do
            pathfound || return assign
        end
    end
    return assign
end
