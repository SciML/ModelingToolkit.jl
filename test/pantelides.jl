using Test
using ModelingToolkit: BiGraph, construct_augmentpath!, init_assign

function augment_path(G)
  assign = init_assign(G)
  pathfound = Ref(false)
  construct_augmentpath!(G, assign, pathfound)
  return assign, pathfound[]
end

# E node using integer from 1 represents the equation mark (the first
# Vector in BiGraph). V node using integer from 11 represents the variables
# (the second vector) Edge in the third part represents their relations

G = BiGraph([1, 2, 3, 4, 5], [11, 12, 13, 14, 15, 16, 17, 18, 19], [1=>15, 2=>16, 3=>17, 4=>18, 3=>19, 4=>19, 5=>0])
#test1 pantelides fig1-b trivial test
@test augment_path(G) == (Dict(11=>nothing,12=>nothing,13=>nothing,14=>nothing,15=>1,16=>2,17=>3,18=>4,19=>nothing), false)

G = BiGraph([1, 2, 3, 4, 5, 6, 7, 8, 9], [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21], [1=>15, 2=>16, 3=>17, 4=>18, 3=>19, 4=>19, 5=>0, 6=>15, 6=>16, 7=>20, 8=>21])
@test augment_path(G) == (Dict(11=>nothing,12=>nothing,13=>nothing,14=>nothing,15=>1,16=>2,17=>3,18=>4,19=>nothing,20=>7,21=>8), false)

G = BiGraph([1, 2, 3, 4, 5, 6, 7, 8, 9], [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21], [1=>15, 2=>16, 3=>17, 4=>18, 3=>19, 4=>19, 5=>0, 7=>20, 8=>21, 7=>17, 8=>18, 9=>20, 9=>21])
#test3 pantelides fig1-d test
@test augment_path(G) == (Dict(11=>nothing,12=>nothing,13=>nothing,14=>nothing,15=>1,16=>2,17=>7,18=>4,19=>3,20=>9,21=>8), true)

G = BiGraph([1, 2, 3, 4], [11, 12, 13, 14, 15, 16], [1=>13, 1=>15, 2=>14, 2=>15, 2=>16, 3=>15])
#test4 pantelides fig2-b test
@test augment_path(G) == (Dict(11=>nothing,12=>nothing,13=>1,14=>2,15=>3,16=>nothing), false)

G = BiGraph([1, 2, 3, 4, 5, 6, 7, 8], [11, 12, 13, 14, 15, 16, 17, 18], [1=>13, 2=>14, 2=>16, 3=>15, 4=>0, 5=>0, 6=>17, 6=>18, 7=>14, 7=>18, 8=>17])
#test5 pantelides fig2-d test
@test augment_path(G) == (Dict(11=>nothing,12=>nothing,13=>1,14=>7,15=>3,16=>2,17=>8,18=>6), true)

G = BiGraph([1, 2, 3, 4], [11, 12, 13, 14], [1=>13, 1=>14, 2=>12, 1=>0])
#test6 pantelides fig3-b trivial test
@test augment_path(G) == (Dict(11=>nothing,12=>2,13=>1,14=>nothing), false)
