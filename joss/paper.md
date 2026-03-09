---
title: '$\text{CO}_2$BatchFill.jl: Efficient multi-layer $\text{CO}_2$ migration modeling using spill-point analysis'
tags:
  - Julia
  - $\text{CO}_2$ storage
  - geological sequestration
  - spill-point analysis
  - reservoir simulation
authors:
  - name: Elling Svee
    orcid: 0009-0008-4225-964X
    corresponding: true
    affiliation: 1
affiliations:
 - name: Norwegian University of Science and Technology, Norway
   index: 1
date: 23 February 2026
bibliography: paper.bib
---

# Summary

When $\text{CO}_2$ is injected into a deep saline aquifer for geological storage, it migrates upward and accumulates beneath low-permeability shale layers in structural traps. If pressure builds sufficiently, $\text{CO}_2$ can breach the shale and enter overlying layers, where it may fill new traps and trigger further migration. `CO2BatchFill.jl` is a Julia package for modeling this type of scenario. Building on the trap detection and filling algorithms in `SurfaceWaterIntegratedModeling.jl` (SWIM), it extends them to handle migration across multiple stacked layers. Simulations run in seconds, which is orders of magnitude less than numerical time-stepping approaches, making the package well suited for workflows that require many model evaluations.

# Statement of need

To model the physical processes behind $\text{CO}_2$ migration, full-physics simulators like OPM Flow [@OPMFlow:2019] use finite-volume methods to solve multi-phase flow equations. Yet, the complexity of the modeling problem and the inherent uncertainty in subsurface properties means different simulators and parameter assumptions can yield drastically different results [@SP11:2025]. One way to handle this uncertainty is to run many simulations with different parameter combinations, but this is often infeasible due to their computational cost. There is therefore a need for interpretable models that can capture the key features of the horizontal and vertical migration, while being fast enough for applications like uncertainty quantification.

# State of the field

`CO2BatchFill.jl` models multi-layer $\text{CO}_2$ migration by combining spill-point analysis with an invasion percolation (IP) approach for vertical migration. Spill-point analysis for $\text{CO}_2$ storage capacity estimation originates from the MATLAB-based MRST-co2lab toolbox [@Nilsen:2015a] and was later reimplemented in Julia as SWIM [@Andersen:2025], though with a focus on surface water flooding rather than subsurface $\text{CO}_2$ modeling. However, both tools only operate on individual layers. Another alternative is the commercial software Permedia [@Permedia:2026], which models multi-layer $\text{CO}_2$ migration using a similar approach. However, as it is not open-source, `CO2BatchFill.jl` is to our knowledge the only open-source package that offers this functionality.

# Principle of invasion percolation

The modeling approach implemented in `CO2BatchFill.jl` is based on IP theory, which assumes that buoyancy and capillary pressure are the governing forces of $\text{CO}_2$ migration in the subsurface [@CallioliInvasion:2025]. This leads to the $\text{CO}_2$ plume filling neighboring regions from lowest to highest capillary entry pressure. Assuming homogeneous properties within each sand layer, the filling order is determined by the topography of the layer, which can be represented through a spillgraph. This spillgraph is a directed graph where nodes represent structural traps and edges represent potential flow paths [@Nilsen:2015a].

Vertical migration occurs when the $\text{CO}_2$ column height in a trap exceeds the capillary entry pressure of the overlying shale. At that point, breaching the shale requires less pressure than filling the next horizontal trap, so the $\text{CO}_2$ migrates vertically [@Carruthers:1998]. The breach location is then treated as a new injection point in the overlying layer, and the same filling algorithm is applied. As the $\text{CO}_2$ drains upward, a fraction is immobilized as residual trapping and left behind the migrating plume. \autoref{fig:ip_system} shows how this process works, and how a graph is used to represent the trap structure.


![Simple example of vertical $\text{CO}_2$ migration with residual trapping. The figure to the left shows the state of a two-layer system before and after migration. A dark blue color is used to represent the $\text{CO}_2$ before drainage, while the light blue color represents the residually trapped $\text{CO}_2$. To the right, the structural traps in the system is represented as nodes in a graph.\label{fig:ip_system}](figures/combined_migration_illustration.svg)

# Software design

<!-- The package is organized around a bottom-to-top layer processing pipeline. The user defines a stack of topography surfaces and reservoir properties, and `analyze_base_surfaces` runs SWIM's `spillanalysis` on each surface to identify structural traps and build a spillgraph — a directed graph encoding the filling order. The main function `fill_layers` then processes layers sequentially from the deepest to the shallowest. Within each layer, traps are filled according to the spillgraph until a $\text{CO}_2$ column exceeds the capillary entry pressure of the overlying shale, at which point the breach is converted into an injection event for the next layer.

A key design choice is to separate the trap-filling algorithm (provided by SWIM) from the multi-layer coupling and leakage logic (implemented in `CO2BatchFill.jl`). This keeps the package lightweight and avoids reimplementing spill-point analysis, while allowing independent updates to either component. Boundary conditions (open or closed) are handled at the topography level by optionally padding surfaces with high-elevation walls before trap analysis.

Visualization is implemented as a package extension that is activated when a Makie backend is loaded. This avoids imposing a heavy graphics dependency on users who only need the simulation functionality. -->


The package is implemented in Julia [@Bezanson:2017], which offers a rich ecosystem of scientific computing. Since many of the algorithms are not easily parallelizable, Julia's JIT compilation are important for achieving high performance.

For setting up a simulation, the user defines a stack of topography surfaces, each with its own set of reservoir properties. For this, the package provides `GenericTopography` and `ReservoirProperties` structs. Having defined the physical properties of the system, the spill-point analysis from SWIM is used on each layer to identify structural traps and build a spillgraph. Defining an injection schedule, the `fill_sequence` function from SWIM can be used to fill traps in each layer. However, to support vertical migration and residual trapping, `CO2BatchFill.jl` extends its functionality in `fill_layer.jl`. Based on the IP approach, the package model the vertical migration by tracking the migration locations and treating them as new injection points in the overlying layer. Constructing an injection schedule that accounts for both the original injection and the migration events, `fill_sequence` can be applied to each layer sequentially from the deepest to the shallowest. 

In addition to the core algorithms, the package also provides utilities for visualizing the $\text{CO}_2$ plumes and their evolution over time. This is implemented as an optional extension that is activated when a Makie backend is loaded. This avoids imposing a heavy graphics dependency on users who only need the simulation functionality.  As a main usecase of the package is uncertainty quantification, an R interface is also provided via the JuliaCall package. This allows users to run simulations and generate figures directly from R.

# Example

We demonstrate `CO2BatchFill.jl` on a synthetic three-layer reservoir, and the full setup is available in `examples/multi_layer_filling.jl`. The reservoir has closed boundaries, meaning $\text{CO}_2$ cannot leak out of the borders of the domain. Each layer is defined on a $100 \times 100$ grid with $10 \mathrm{m}$ cell spacing, and sampled from a Gaussian random field with a Matérn covariance function. The two bottom shale layers have a capillary entry pressure of $25 \text{kPa}$, while the top layer is impermeable caprock. The sand layers have a residual trapping of $40 \%$, while the rest will drain out over $5$ years.

$\text{CO}_2$ is injected at a constant rate of $80 \text{m}^3/\text{year}$ over $10$ years into the center of the bottom layer. It takes approximately $5.5$ years for the first breach to occur. \autoref{fig:co2_L1_plume} shows the $\text{CO}_2$ plume in the first layer at three steps up to the time of migration. As more and more $\text{CO}_2$ is injected, the plume grows and fills traps until it eventually reaches the capillary entry pressure of the overlying shale. \autoref{fig:co2_plume} shows the plumes in the three layers after $15$ years. At this point, the $\text{CO}_2$ has migrated all the way to the top layer. There, the impermeable caprock prevents further migration, resulting in a buildup of $\text{CO}_2$. 

![$\text{CO}_2$ plume in storage layer $1$ for three time steps up to the time of migration. The injection location is marked by a cross.\label{fig:co2_L1_plume}](figures/single_layer_snapshots_L1.svg){width=100%}

![$\text{CO}_2$ plume extents after $15$ years. The injection location is marked by a cross, and the leakage locations are indicated by triangles.\label{fig:co2_plume}](figures/multi_layer_co2_final.svg){width=100%}


To investigate the impact of how different capillary entry pressures affect the migration, we can run an ensemble of simulations with different values. This is often done in practice to account for uncertainty in subsurface properties. We run an ensemble of $100$ simulations with entry pressures varying between $20$ and $30 \text{kPa}$. The total runtime of the full ensemble is only $9$ seconds, and the implementation is available in `examples/multi_layer_ensemble.jl`.

\autoref{fig:ensemble_timeseries} shows the stored and drained $\text{CO}_2$ volumes in each layer over time. The lines show the mean across the ensemble, while the shaded areas indicate the standard deviation. Observe an initial increase in stored $\text{CO}_2$ as traps fill, followed by a decrease once migration begins. Due to residual trapping, the total $\text{CO}_2$ mass does not drop to zero, but stabilizes at a certain level. The variability across the ensemble is relatively large, which highlights the importance of uncertainty quantification in $\text{CO}_2$ storage modeling.


![Ensemble of $100$ simulations with varying capillary entry pressures. The lines show the mean of the total stored and drained $\text{CO}_2$ across the ensemble, while the shaded areas indicate the standard deviation. \label{fig:ensemble_timeseries}](figures/ensemble_timeseries.svg){width=100%}

# Research impact statement

The package was originally developed as a surrogate model for uncertainty quantification in $\text{CO}_2$ storage. Ongoing research applies it for analyzing the Sleipner storage site, where it is used to explore how different parameter assumptions affect storage capacity estimates. It is also being used for Bayesian optimization of $\text{CO}_2$ injection strategies, where its speed allows for efficient optimization loops that would be infeasible with full-physics simulators.

# AI usage disclosure

The main architectural design and modeling approach was created by the authors without AI assistance.

Claude Code Sonnet $4.5$ and Opus $4.6$ were used for code generation during development. The authors used these tools to generate code snippets, documentation, and to assist with debugging and refactoring. This was done iteratively and in conjunction with manual coding, with the authors providing detailed prompts for specific tasks. All generated content was reviewed, validated and edited by the authors. GitHub Copilot helped with autocompletion. 

The paper was written by the authors without AI assistance, though the authors used ChatGPT $5.3$ to check and improve grammar.

# Acknowledgements

We thank Odd A. Andersen (SINTEF) for developing SWIM, which provides the spill-point analysis algorithms that `CO2BatchFill.jl` builds upon. Also thanks to Philip Ringrose for his insights on $\text{CO}_2$ storage, which helped shape the modeling approach and features of the package.


# References
