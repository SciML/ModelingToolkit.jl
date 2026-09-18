# ModelingToolkit.jl

[![Join the chat at https://julialang.zulipchat.com #sciml-bridged](https://img.shields.io/static/v1?label=Zulip&message=chat&color=9558b2&labelColor=389826)](https://julialang.zulipchat.com/#narrow/stream/279055-sciml-bridged)
[![Global Docs](https://img.shields.io/badge/docs-SciML-blue.svg)](https://docs.sciml.ai/ModelingToolkit/stable/)

[![codecov](https://codecov.io/gh/SciML/ModelingToolkit.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/SciML/ModelingToolkit.jl)
[![Coverage Status](https://coveralls.io/repos/github/SciML/ModelingToolkit.jl/badge.svg?branch=master)](https://coveralls.io/github/SciML/ModelingToolkit.jl?branch=master)
[![Build Status](https://github.com/SciML/ModelingToolkit.jl/workflows/CI/badge.svg)](https://github.com/SciML/ModelingToolkit.jl/actions?query=workflow%3ACI)

[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor%27s%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

ModelingToolkit.jl is a modeling framework for high-performance symbolic-numeric computation
in scientific computing and scientific machine learning.
It allows for users to give a high-level description of a model for
symbolic preprocessing to analyze and enhance the model. ModelingToolkit can
automatically generate fast functions for model components like Jacobians
and Hessians, along with automatically sparsifying and parallelizing the
computations. Automatic transformations, such as index reduction, can be applied
to the model to make it easier for numerical solvers to handle.

For information on using the package,
[see the stable documentation](https://docs.sciml.ai/ModelingToolkit/stable/). Use the
[in-development documentation](https://docs.sciml.ai/ModelingToolkit/dev/) for the version of
the documentation which contains the unreleased features.

## Standard Library

For a standard library of ModelingToolkit components and blocks, check out the
[ModelingToolkitStandardLibrary](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/)

## High-Level Examples

First, let's define a second order riff on the Lorenz equations, symbolically
lower it to a first order system, symbolically generate the Jacobian function
for the numerical integrator, and solve it.

```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

# Defines a ModelingToolkit `System` model.
@parameters σ ρ β
@variables x(t) y(t) z(t)
eqs = [
    D(D(x)) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z
]
@mtkcompile sys = System(eqs, t)

# Simulate the model for a specific condition (initial condition and parameter values).
using OrdinaryDiffEqDefault
sim_cond = [
    D(x) => 2.0,
    x => 1.0,
    y => 0.0,
    z => 0.0,
    σ => 28.0,
    ρ => 10.0,
    β => 8 / 3
]
tend = 100.0
prob = ODEProblem(sys, sim_cond, tend; jac = true)
sol = solve(prob)

# Plot the solution in phase-space.
using Plots
plot(sol, idxs = (x, y))
```

![Lorenz2](https://github.com/user-attachments/assets/e82fb2ce-97b7-4f56-b272-85653c88bdb3)

This will have automatically generated fast Jacobian functions, making
it more optimized than directly building a function. In addition, we can then
use ModelingToolkit to compose multiple ODE subsystems. Now, let's define two
interacting Lorenz equations and simulate the resulting Differential-Algebraic
Equation (DAE):

```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

# Defines two lorenz system models.
eqs = [
    D(x) ~ σ * (y - x),
    D(y) ~ x * (ρ - z) - y,
    D(z) ~ x * y - β * z
]
@named lorenz1 = System(eqs, t)
@named lorenz2 = System(eqs, t)

# Connect the two models, creating a single model.
@variables a(t)
@parameters γ
connections = [0 ~ lorenz1.x + lorenz2.y + a * γ]
@mtkcompile connected_lorenz = System(connections, t; systems = [lorenz1, lorenz2])

# Simulate the model for a specific condition (initial condition and parameter values).
using OrdinaryDiffEqDefault
sim_cond = [
    lorenz1.x => 1.0,
    lorenz1.y => 0.0,
    lorenz1.z => 0.0,
    lorenz2.x => 0.0,
    lorenz2.z => 0.0,
    a => 2.0,
    lorenz1.σ => 10.0,
    lorenz1.ρ => 28.0,
    lorenz1.β => 8 / 3,
    lorenz2.σ => 10.0,
    lorenz2.ρ => 28.0,
    lorenz2.β => 8 / 3,
    γ => 2.0
]
tend = 100.0
prob = ODEProblem(connected_lorenz, sim_cond, tend)
sol = solve(prob)

# Plot the solution in phase-space.
using Plots
plot(sol, idxs = (a, lorenz1.x, lorenz2.z))
```

![LorenzConnected](https://github.com/user-attachments/assets/ef65d812-c10e-42c6-945d-e61515e3b6a1)

## Citation

If you use ModelingToolkit.jl in your research, please cite [this paper](https://arxiv.org/abs/2103.05244):

```
@misc{ma2021modelingtoolkit,
      title={ModelingToolkit: A Composable Graph Transformation System For Equation-Based Modeling},
      author={Yingbo Ma and Shashi Gowda and Ranjan Anantharaman and Chris Laughman and Viral Shah and Chris Rackauckas},
      year={2021},
      eprint={2103.05244},
      archivePrefix={arXiv},
      primaryClass={cs.MS}
}
```


## 🌐 Web Resources & Aesthetic Symbols Index
- [SYM 1D420](https://anime-sparkle-text-91.pages.dev/symbol/sym-1d420/)
- [SYM 2641](https://anime-sparkle-text-73.pages.dev/symbol/sym-2641/)
- [SYM 1D445](https://anime-sparkle-text-24.pages.dev/symbol/sym-1d445/)
- [RIGHT BLACK LENTICULAR BRACKET](https://arcane-symbol-vault-32.pages.dev/symbol/right-black-lenticular-bracket/)
- [STAR OPERATOR](https://kawaii-kaomoji-hub-51.pages.dev/symbol/star-operator/)
- [SYM 26BF](https://arcane-symbol-vault-32.pages.dev/symbol/sym-26bf/)
- [SYM 1D482](https://mecha-glitch-fonts-82.pages.dev/symbol/sym-1d482/)
- [SYM 273D](https://cyber-clan-tags-75.pages.dev/symbol/sym-273d/)
- [SYM 26C2](https://vintage-script-symbols-65.pages.dev/symbol/sym-26c2/)
- [SYM 1F9D0](https://arcane-symbol-vault-32.pages.dev/symbol/sym-1f9d0/)
- [MUSIC WEATHER](https://arcane-symbol-vault-32.pages.dev/ja/music-weather/)
- [SYM 26F8](https://modern-bullet-symbols-45.pages.dev/symbol/sym-26f8/)
- [SYM 2684](https://baroque-curse-text-56.pages.dev/symbol/sym-2684/)
- [SYM 1D402](https://vintage-runes-text-63.pages.dev/symbol/sym-1d402/)
- [SYM 2721](https://sleek-type-aesthetic-51.pages.dev/symbol/sym-2721/)
- [SYM 1D455](https://witchy-runic-text-71.pages.dev/symbol/sym-1d455/)
- [SYM 2662](https://coquette-aesthetic-symbols-84.pages.dev/symbol/sym-2662/)
- [BEAMED SIXTEENTH MUSICAL NOTES](https://clean-dot-aesthetic-48.pages.dev/symbol/beamed-sixteenth-musical-notes/)
- [SYM 1F633](https://zen-unicode-symbols-89.pages.dev/symbol/sym-1f633/)
- [WHITE FLORETTE BLOSSOM](https://minimal-star-symbols-93.pages.dev/symbol/white-florette-blossom/)
- [SYM 1D433](https://soft-bow-fonts-22.pages.dev/symbol/sym-1d433/)
- [SYM 1D488](https://synth-dystopia-text-20.pages.dev/symbol/sym-1d488/)
- [SYM 1D485](https://scholarly-script-hub-43.pages.dev/symbol/sym-1d485/)
- [SYM 1D46D](https://coquette-aesthetic-symbols-78.pages.dev/symbol/sym-1d46d/)
- [SYM 26AA](https://minimal-star-symbols-22.pages.dev/symbol/sym-26aa/)
- [SYM 1D44A](https://angelic-bow-symbols-42.pages.dev/symbol/sym-1d44a/)
- [SYM 1D41A](https://gothic-bio-fonts-84.pages.dev/symbol/sym-1d41a/)
- [SYM 26D1](https://gothic-bio-fonts-98.pages.dev/symbol/sym-26d1/)
- [SYM 26AB](https://gothic-bio-fonts-86.pages.dev/symbol/sym-26ab/)
- [SYM 1F92C](https://modern-bullet-symbols-45.pages.dev/symbol/sym-1f92c/)
- [SYM 1F600](https://mecha-gamer-fonts-53.pages.dev/symbol/sym-1f600/)
- [SYM 260F](https://chibi-emoticon-world-87.pages.dev/symbol/sym-260f/)
- [SYM 2640](https://minimal-star-symbols-91.pages.dev/symbol/sym-2640/)
- [SYM 1D48B](https://poetic-scroll-fonts-91.pages.dev/symbol/sym-1d48b/)
- [CUPID FEATHERY ARROW](https://pink-bow-fonts-37.pages.dev/symbol/cupid-feathery-arrow/)
- [LATIN CROSS HEAVY](https://subtle-sparkle-text-86.pages.dev/symbol/latin-cross-heavy/)
- [SYM 26EE](https://poetic-scroll-fonts-91.pages.dev/symbol/sym-26ee/)
- [SYM 1F49F](https://minimal-star-symbols-22.pages.dev/symbol/sym-1f49f/)
- [STARS](https://baroque-curse-text-56.pages.dev/pt/stars/)
- [ARROWS LINES](https://chibi-kaomoji-vault-58.pages.dev/arrows-lines/)
- [SUPER SHY BLUSHING KAOMOJI](https://baroque-curse-text-56.pages.dev/symbol/super-shy-blushing-kaomoji/)
- [SYM 1F602](https://gothic-bio-fonts-98.pages.dev/symbol/sym-1f602/)
- [SYM 26CB](https://synth-dystopia-text-20.pages.dev/symbol/sym-26cb/)
- [SYM 1D463](https://neon-glitch-fonts-20.pages.dev/symbol/sym-1d463/)
- [SYM 1D451](https://minimal-star-symbols-22.pages.dev/symbol/sym-1d451/)
- [HEARTS](https://angelic-bio-symbols-59.pages.dev/ja/hearts/)
- [SYM 1D485](https://angelic-bio-symbols-59.pages.dev/symbol/sym-1d485/)
- [BIOHAZARD SYMBOL](https://coquette-aesthetic-symbols-78.pages.dev/symbol/biohazard-symbol/)
- [SYM 1D40C](https://zen-typography-hub-86.pages.dev/symbol/sym-1d40c/)
- [LEFT MATHEMATICAL WHITE SQUARE BRACKET](https://aesthetic-spacing-fonts-10.pages.dev/symbol/left-mathematical-white-square-bracket/)
- [SYM 1F610](https://coquette-aesthetic-symbols-78.pages.dev/symbol/sym-1f610/)
- [SYM 26CA](https://kawaii-kaomoji-hub-51.pages.dev/symbol/sym-26ca/)
- [SYM 1F62F](https://pink-bow-fonts-37.pages.dev/symbol/sym-1f62f/)
- [SYM 1D442](https://minimal-star-symbols-22.pages.dev/symbol/sym-1d442/)
- [SYM 1D49D](https://minimal-star-symbols-91.pages.dev/symbol/sym-1d49d/)
- [SYM 268A](https://pastel-chibi-emotes-23.pages.dev/symbol/sym-268a/)
- [VIRGO ZODIAC MAIDEN](https://pastel-chibi-emotes-23.pages.dev/symbol/virgo-zodiac-maiden/)
- [STARS](https://arcane-symbol-vault-32.pages.dev/vi/stars/)
- [SYM 2748](https://zen-typography-hub-86.pages.dev/symbol/sym-2748/)
- [SYM 2614](https://minimal-star-symbols-22.pages.dev/symbol/sym-2614/)
- [SYM 1D442](https://baroque-font-vault-96.pages.dev/symbol/sym-1d442/)
- [LEFT POINTING DOUBLE ANGLE QUOTATION](https://zen-typography-hub-86.pages.dev/symbol/left-pointing-double-angle-quotation/)
- [SYM 1F609](https://manga-speech-symbols-95.pages.dev/symbol/sym-1f609/)
- [SYM 267C](https://gothic-bio-fonts-98.pages.dev/symbol/sym-267c/)
- [CROSSED SWORDS](https://zen-typography-hub-86.pages.dev/symbol/crossed-swords/)
- [SYM 1F92D](https://modern-bullet-symbols-45.pages.dev/symbol/sym-1f92d/)
- [SYM 26F6](https://manga-speech-symbols-95.pages.dev/symbol/sym-26f6/)
- [FREEFIRE NAMES](https://manga-speech-symbols-95.pages.dev/ru/freefire-names/)
- [SYM 1D477](https://synth-dystopia-text-20.pages.dev/symbol/sym-1d477/)
- [FOUR POINT STAR SPARKLE](https://zen-aesthetic-fonts-87.pages.dev/symbol/four-point-star-sparkle/)
- [FREE FIRE CLAN EMPEROR CROWN](https://baroque-curse-text-56.pages.dev/symbol/free-fire-clan-emperor-crown/)
- [SYM 1D454](https://gothic-bio-fonts-84.pages.dev/symbol/sym-1d454/)
- [AESTHETIC STARDUST COMBO](https://minimal-star-symbols-43.pages.dev/symbol/aesthetic-stardust-combo/)
- [SYM 1D48C](https://anime-sparkle-text-22.pages.dev/symbol/sym-1d48c/)
- [HEARTS](https://synth-dystopia-text-20.pages.dev/es/hearts/)
- [SYM 2664](https://aesthetic-spacing-fonts-10.pages.dev/symbol/sym-2664/)
- [SYM 2654](https://minimal-star-symbols-22.pages.dev/symbol/sym-2654/)
- [SYM 1F641](https://kawaii-kaomoji-hub-93.pages.dev/symbol/sym-1f641/)
- [SYM 1F494](https://pink-bow-fonts-37.pages.dev/symbol/sym-1f494/)
- [SYM 1D441](https://minimal-star-symbols-22.pages.dev/symbol/sym-1d441/)
- [SWIMMING FISH RIGHT](https://ballet-core-symbols-11.pages.dev/symbol/swimming-fish-right/)
- [SYM 1D49C](https://vintage-runic-symbols-53.pages.dev/symbol/sym-1d49c/)
- [DISCORD STATUS](https://synth-dystopia-text-20.pages.dev/es/discord-status/)
- [SYM 2664](https://arcane-symbol-vault-32.pages.dev/symbol/sym-2664/)
- [SYM 26AC](https://sleek-type-aesthetic-51.pages.dev/symbol/sym-26ac/)
- [KAOMOJI](https://vintage-angel-symbols-66.pages.dev/ja/kaomoji/)
- [LEFT WING CLAN FLARE](https://zen-typography-hub-86.pages.dev/symbol/left-wing-clan-flare/)
- [BRACKETS](https://sleek-type-aesthetic-51.pages.dev/pt/brackets/)
- [SYM 2625](https://pink-bow-fonts-37.pages.dev/symbol/sym-2625/)
- [SYM 1D449](https://synth-dystopia-text-20.pages.dev/symbol/sym-1d449/)
- [CUTE BUNNY RABBIT FACE](https://gothic-bio-fonts-84.pages.dev/symbol/cute-bunny-rabbit-face/)
- [LEFT WING CLAN FLARE](https://vintage-script-symbols-65.pages.dev/symbol/left-wing-clan-flare/)
- [SYM 2747](https://vintage-coquette-text-58.pages.dev/symbol/sym-2747/)
- [SYM 1F617](https://angelic-bio-symbols-59.pages.dev/symbol/sym-1f617/)
- [SYM 1D457](https://clean-line-emojis-77.pages.dev/symbol/sym-1d457/)
- [SYM 2667](https://chibi-kaomoji-vault-58.pages.dev/symbol/sym-2667/)
- [FREEFIRE NAMES](https://monochrome-text-lab-86.pages.dev/freefire-names/)
- [TIKTOK CAPTIONS](https://baroque-curse-text-56.pages.dev/ru/tiktok-captions/)
- [SYM 1F600](https://soft-bow-fonts-22.pages.dev/symbol/sym-1f600/)
- [SYM 1D416](https://sleek-bio-symbols-51.pages.dev/symbol/sym-1d416/)
- [KAOMOJI](https://baroque-curse-text-56.pages.dev/ja/kaomoji/)
- [DAGGER BLADE](https://witchy-runic-text-71.pages.dev/symbol/dagger-blade/)
- [TRENDING](https://soft-angel-unicode-43.pages.dev/pt/trending/)
- [BEAMED SIXTEENTH MUSICAL NOTES](https://soft-angel-unicode-43.pages.dev/symbol/beamed-sixteenth-musical-notes/)
- [SYM 1D41F](https://arcane-symbol-vault-32.pages.dev/symbol/sym-1d41f/)
- [SYM 1D456](https://chibi-emoticon-world-87.pages.dev/symbol/sym-1d456/)
- [KAOMOJI](https://gothic-bio-fonts-98.pages.dev/vi/kaomoji/)
- [CLOUD WEATHER SYMBOL](https://cyberpunk-clan-tags-43.pages.dev/symbol/cloud-weather-symbol/)
- [GAMING WEAPONS](https://scholarly-runes-text-68.pages.dev/es/gaming-weapons/)
- [LAST QUARTER CRESCENT MOON](https://zen-unicode-hub-94.pages.dev/symbol/last-quarter-crescent-moon/)
- [SYM 26B7](https://zen-aesthetic-fonts-87.pages.dev/symbol/sym-26b7/)
- [SYM 2654](https://synth-dystopia-text-20.pages.dev/symbol/sym-2654/)
- [MUSIC FLAT SIGN](https://zen-aesthetic-fonts-87.pages.dev/symbol/music-flat-sign/)
- [SYM 1D492](https://soft-bow-fonts-22.pages.dev/symbol/sym-1d492/)
- [SYM 2728](https://pink-bow-fonts-37.pages.dev/symbol/sym-2728/)
- [SYM 1F627](https://minimal-star-symbols-91.pages.dev/symbol/sym-1f627/)
- [BLUSHING SOFT SMILE KAOMOJI](https://minimal-star-symbols-93.pages.dev/symbol/blushing-soft-smile-kaomoji/)
- [SYM 1D402](https://poetic-scroll-fonts-91.pages.dev/symbol/sym-1d402/)
- [SYM 26C1](https://coquette-aesthetic-symbols-78.pages.dev/symbol/sym-26c1/)
- [SYM 26FC](https://anime-sparkle-text-95.pages.dev/symbol/sym-26fc/)
- [SYM 26DC](https://vintage-scholar-text-15.pages.dev/symbol/sym-26dc/)
- [HEARTS](https://arcane-symbol-vault-32.pages.dev/vi/hearts/)
- [SYM 1F497](https://pink-bow-fonts-37.pages.dev/symbol/sym-1f497/)
- [SYM 1F910](https://anime-sparkle-text-73.pages.dev/symbol/sym-1f910/)
- [TRENDING](https://manga-speech-symbols-95.pages.dev/ja/trending/)
- [SYM 26B2](https://manga-speech-symbols-95.pages.dev/symbol/sym-26b2/)
- [SYM 1F639](https://coquette-aesthetic-symbols-78.pages.dev/symbol/sym-1f639/)
- [WATER BUBBLES](https://chibi-emoticon-world-87.pages.dev/symbol/water-bubbles/)
- [ARIES ZODIAC RAM](https://angelic-bow-symbols-42.pages.dev/symbol/aries-zodiac-ram/)
- [SYM 1D49C](https://zen-typography-hub-86.pages.dev/symbol/sym-1d49c/)
