using Test
using ModelingToolkit: BiGraph, augment_path
# E node using integer from 1 represents the equation mark (the first
# Vector in BiGraph). V node using integer from 11 represents the variables
# (the second vector) Edge in the third part represents their relations

G=BiGraph(Vector([1,2,3,4,5]),Vector([11,12,13,14,15,16,17,18,19]),[[1,15],[2,16],[3,17],[4,18],[3,19],[4,19],[5,0]])

#test1 pantelides fig1-b trivial test
@test augment_path(G) == Dict(11=>0,12=>0,13=>0,14=>0,15=>1,16=>2,17=>3,18=>4,19=>0)

G=BiGraph(Vector([1,2,3,4,5,6,7,8,9]),Vector([11,12,13,14,15,16,17,18,19,20,21]),[[1,15],[2,16],[3,17],[4,18],[3,19],[4,19],[5,0],[6,15],[6,16],[7,20],[8,21]])
@test augment_path(G) == Dict(11=>0,12=>0,13=>0,14=>0,15=>1,16=>2,17=>3,18=>4,19=>0,20=>7,21=>8)

G=BiGraph(Vector([1,2,3,4,5,6,7,8,9]),Vector([11,12,13,14,15,16,17,18,19,20,21]),[[1,15],[2,16],[3,17],[4,18],[3,19],[4,19],[5,0],[7,20],[8,21],[7,17],[8,18],[9,20],[9,21]])
#test3 pantelides fig1-d  test
@test augment_path(G) == Dict(11=>0,12=>0,13=>0,14=>0,15=>1,16=>2,17=>7,18=>4,19=>3,20=>9,21=>8)

G=BiGraph(Vector([1,2,3,4]),Vector([11,12,13,14,15,16]),[[1,13],[1,15],[2,14],[2,15],[2,16],[3,15]])#[1,1],[1,3],[2,2],[2,4],[2,5],[2,6],[3,1],[3,2],[3,5],[4,1]
#test4 pantelides fig2-b  test
@test augment_path(G) == Dict(11=>0,12=>0,13=>1,14=>2,15=>3,16=>0)

G=BiGraph(Vector([1,2,3,4,5,6,7,8]),Vector([11,12,13,14,15,16,17,18]),[[1,13],[2,14],[2,16],[3,15],[4,0],[5,0],[6,17],[6,18],[7,14],[7,18],[8,17]])#[1,1],[1,3],[2,2],[2,4],[2,5],[2,6],[3,1],[3,2],[3,5],[4,1]
#test5 pantelides fig2-d  test
@test augment_path(G) == Dict(11=>0,12=>0,13=>1,14=>7,15=>3,16=>2,17=>8,18=>6)

G=BiGraph(Vector([1,2,3,4]),Vector([11,12,13,14]),[[1,13],[1,14],[2,12],[1,0]])
#test6 pantelides fig3-b trivial test
@test augment_path(G) == Dict(11=>0,12=>2,13=>1,14=>0)
