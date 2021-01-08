using SpecialFunctions
import Base.Broadcast


const linearity_known_1 = IdDict{Function,Bool}()
const linearity_known_2 = IdDict{Function,Bool}()

const linearity_map_1 = IdDict{Function, Bool}()
const linearity_map_2 = IdDict{Function, Tuple{Bool, Bool, Bool}}()

# 1-arg

const monadic_linear = [deg2rad, +, rad2deg, transpose, -, conj]

const monadic_nonlinear = [asind, log1p, acsch, erfc, digamma, acos, asec, acosh, airybiprime, acsc, cscd, log, tand, log10, csch, asinh, airyai, abs2, gamma, lgamma, erfcx, bessely0, cosh, sin, cos, atan, cospi, cbrt, acosd, bessely1, acoth, erfcinv, erf, dawson, inv, acotd, airyaiprime, erfinv, trigamma, asecd, besselj1, exp, acot, sqrt, sind, sinpi, asech, log2, tan, invdigamma, airybi, exp10, sech, erfi, coth, asin, cotd, cosd, sinh, abs, besselj0, csc, tanh, secd, atand, sec, acscd, cot, exp2, expm1, atanh]

# We store 3 bools even for 1-arg functions for type stability
const three_trues = (true, true, true)
for f in monadic_linear
    linearity_known_1[f] = true
    linearity_map_1[f] = true
end

for f in monadic_nonlinear
    linearity_known_1[f] = true
    linearity_map_1[f] = false
end

# 2-arg
for f in [+, rem2pi, -, >, isless, <, isequal, max, min, convert, <=, >=]
    linearity_known_2[f] = true
    linearity_map_2[f] = (true, true, true)
end

for f in [*]
    linearity_known_2[f] = true
    linearity_map_2[f] = (true, true, false)
end

for f in [/]
    linearity_known_2[f] = true
    linearity_map_2[f] = (true, false, false)
end
for f in [\]
    linearity_known_2[f] = true
    linearity_map_2[f] = (false, true, false)
end

for f in [hypot, atan, mod, rem, lbeta, ^, beta]
    linearity_known_2[f] = true
    linearity_map_2[f] = (false, false, false)
end

haslinearity_1(@nospecialize(f)) = get(linearity_known_1, f, false)
haslinearity_2(@nospecialize(f)) = get(linearity_known_2, f, false)

linearity_1(@nospecialize(f)) = linearity_map_1[f]
linearity_2(@nospecialize(f)) = linearity_map_2[f]

# TermCombination datastructure

struct TermCombination
    terms::Set{Dict{Int, Int}} # idx => pow
end

@eval Base.one(::Type{TermCombination}) = $(TermCombination(Set([Dict{Int,Int}()])))
@eval Base.zero(::Type{TermCombination}) = $(TermCombination(Set{Dict{Int,Int}}()))

#=
function Base.:(==)(comb1::TermCombination, comb2::TermCombination)
    comb1.terms == comb2.terms && return true

    n1 = reduce(max, (k for (k,_) in Iterators.flatten(comb1.terms)), init=0)
    n2 = reduce(max, (k for (k,_) in Iterators.flatten(comb2.terms)), init=0)
    n = max(n1, n2)

    _sparse(comb1, n) == _sparse(comb2, n)
end
=#

# to make Mul and Add work
Base.:*(::Number, comb::TermCombination) = comb
function Base.:^(comb::TermCombination, ::Number)
    isone(comb) && return comb
    iszero(comb) && return _scalar
    return comb *  comb
end

function Base.:+(comb1::TermCombination, comb2::TermCombination)
    if isone(comb1) && !iszero(comb2)
        return comb2
    elseif isone(comb2) && !iszero(comb1)
        return comb1
    elseif comb1 === comb2
        return comb1
    end
    TermCombination(union(comb1.terms, comb2.terms))
end

Base.:+(comb1::TermCombination) = comb1

function _merge(dict1, dict2)
    d = copy(dict1)
    for (k, v) in dict2
        d[k] = min(2, get(dict1, k, 0) + v)
    end
    d
end

function Base.:*(comb1::TermCombination, comb2::TermCombination)
    if isone(comb1)
        return comb2
    elseif isone(comb2)
        return comb1
    elseif comb1 === comb2 # squaring optimization
        terms = comb1.terms
        # turns out it's enough to track
        # a^2*b^2
        # and a^2 + b^2 + ab
        # have the same hessian sparsity
        t = Dict(k=>2 for (k,_) in
                 Iterators.flatten(terms))
        TermCombination(Set([t]))
        #=
        # square each term
        t1 = [Dict(k=>2 for (k,_) in dict)
              for dict in comb1.terms]
        # multiply each term
        t2 = Dict{Int,Int}[]
        for i in 1:length(terms)
            for j in i+1:length(terms)
                push!(t2, _merge(terms[i], terms[j]))
            end
        end
        TermCombination(union(t1, t2))
        =#
    else
        Set([_merge(dict1, dict2)
             for dict1 in comb1.terms,
             dict2 in comb2.terms]) |> TermCombination
    end
end
Base.:*(comb1::TermCombination) = comb1
Base.iszero(c::TermCombination) = isempty(c.terms)
Base.isone(c::TermCombination) = all(isempty, c.terms)

function _sparse(t::TermCombination, n)
    I = Int[]
    J = Int[]
    for dict in t.terms
        kv = collect(pairs(dict))
        for i in 1:length(kv)
            k, v = kv[i]
            if v>=2
                push!(I, k)
                push!(J, k)
            end
            for j in i+1:length(kv)
                if v >= 1 && kv[j][2] >= 1
                    push!(I, k)
                    push!(J, kv[j][1])
                end
            end
        end
    end
    s1 = sparse(I,J,fill!(BitVector(undef, length(I)), true),n,n)
    s1 .| s1'
end

# 1-arg functions
combine_terms_1(lin, term) = lin ? term : term * term

# 2-arg functions
function combine_terms_2(linearity, term1, term2)

    linear11, linear22, linear12 = linearity
    term = zero(TermCombination)
    if linear11
        if !linear12
            term += term1
        end
    else
        term += term1 * term1
    end

    if linear22
        if !linear12
            term += term2
        end
    else
        term += term2 * term2
    end

    if linear12
        term += term1 + term2
    else
        term += term1 * term2
    end
    term
end
