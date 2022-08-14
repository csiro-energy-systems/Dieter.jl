# Model formulation

```@contents
Pages = ["optmodel.md"]
Depth = 3
```

## The optimization model 

The Dieter.jl package contains Julia functions to generate [JuMP](https://github.com/jump-dev/JuMP.jl) linear programming model instances. Each model has a system cost-minimising objective, with  energy and power variables representing technology capacities and operational states. 

While there are degrees of flexibility around the final problem formulation, there are some core structures that remain invariant across models.

### Index sets

Each optimization model is structured around some key indexing sets related to
- Temporal
- Spatial
- Transmission
- Technological


#### Temporal

The time indexing of DIETER is typically over 8760 hours of a single full year. The start and end day of the year itself is not assumed.

For example, in the code we have something akin to
```julia
hoursInYear = 8760
Timestep = 1
Hours = 1:Timestep:hoursInYear
# [1, 2, ..., 8759, 8760]
```

#### Spatial

The model is represented spatially as a network graph, where generation and demand activity is located at nodes (or vertices), while connections between nodes are represented as directed edges (or lines). A typical node in the model represents a region of a system in which we aggregate the operations and capacities of technologies. We may assign to a node a demand requirement, for electricity or another commodity, that must be satisfied by energy production.

For example, we might define a set of nodes as a simple array of string labels:
```julia
Nodes = ["R1", "R2", "R3"]
```

We use both the terms "region" and "zone" to refer to spatial entities, often using "region" to refer to the largest spatial grouping in the model, and a "zone" to mean a sub-region.
Specific examples are transmission zones (abbreviated `TxZones`) and renewable energy zones (abbreviated `REZones`, used for aggregating variable renewable energy resources).

#### Transmission

Transmission lines transfer electrical power from supply to demand. As mentioned, regions or zones are connected by lines, represented as ordered pairs of spatial nodes. We refer to these lines as `Arcs` in the models: an arc may represent what in reality is a transmission link comprised of one or more physical transmission lines.

A basic line representation is a pair consisting of two spatial regions. We collect all lines in the set `Arcs`

For example:
```julia
Arcs = [("R1","R2"), ("R2", "R1"), ("R2,R3"), ("R2", "R3")]
```
Here the network topology is a linear graph (R1 ⇆ R2 ⇆ R3) where both directions of a connection between regions are explicitly included in the set.

#### Technological

Technologies are sub-categorised as follows:
- Conventional (fossil fuel-based) electrical generators
- Variable renewable electrical generators
- Energy storage devices
- Hydrogen production devices

The set `Techs` contains variable and dispatchable electricity generation types, for example
```julia
Techs = ["Wind", "SolarPV", "Hydro", "CCGT", "OCGT"]
```
and is distinguished from the set `Storages` of energy storage technology types, for example
```julia
Storages = ["LIonBattery", "PumpedHydro", "CAES", "RedoxBattery"]
```
since the model constraints differ for these two categories.

#### Mapping between index sets

Any model must construct relations between basic components in different categories.

For instance, it is important to consider that not all technology types exist in every spatial region. It is therefore necessary to have a method of mapping technologies to regions. A direct method is to maintain an index array of pairs consisting of a node matched with a technology:
```julia
Nodes_Techs = [ ("R1","Hydro"), ("R1","CCGT"),("R1","SolarPV"), 
                ("R2","SolarPV"), ("R2","Wind"),
                ("R3","Hydro"), ("R3","CCGT"), ("R3","SolarPV"), ("R3","Wind")]
``` 
while for storage technologies the correspondence between location and available storage technologies is similarly defined through the array of pairs
```julia
Nodes_Storages = [ ("R1","LIonBattery"),
                   ("R3","LIonBattery"), ("R3","PumpedHydro")]
```

Numerous other similar subsets are used to define sets of constraints on particular combinations of technology, spatial location and so on. For example, `Nodes_Dispatch` is the subset of `Nodes_Techs` filtered to have only dispatchable technologies in the second component.

### Variables

The key decision variables of the model are shown here as

 `Julia variable name` - *Short descriptor* - ``\TeX`` notation
 
 The notation `v[I,J]` indicates that the variable `v` is jointly indexed over `I`, `J`, etc.

#### Generation operational variables
- `G[Nodes_Techs, Hours]` - *Generation level* - ``G_{(n,t),h}``
    - Detail: Generation level - all generation technologies
    - Bounds: lower_bound=0
    - Units: MWh per time-interval; 
- `G_UP[Nodes_Dispatch, Hours]` - *Generation upshift* - ``G^{\uparrow}_{(n,t),h}``
    - Detail: Generation level change upwards
    - Bound: lower_bound=0
    - Units: MWh per time-interval; 
- `G_DO[Nodes_Dispatch, Hours]` - *Generation downshift* - ``G^{\downarrow}_{(n,t),h}``
    - Detail: Generation level change downwards
    - Bound: lower_bound=0
    - Units: MWh per time-interval; 
- `G_REZ[REZones,Hours]` - *Generation renewable* - ``G^{\text{REZ}}_{z,h}``
    - Detail: Generation level - renewable energy zone tech. & storage
    - Bound: lower_bound=0
    - Units: MWh per time-interval; 
    - __Note__: This is a 'book-keeping' variable for each renewable energy zone (REZ, `REZones`) that aggregates operational variables at each hour.
- `G_TxZ[TxZones,Hours]` - *Generation transmission zones* - ``G^{\text{TxZ}}_{z,h}``
    - Detail: Generation level - transmission zone tech. & storage
    - Bound: lower_bound=0
    - Units: MWh per time-interval; 
    - __Note__: This is a 'book-keeping' variable for each transmission zone (`TxZones`) that aggregates operational variables in each hour.

#### Storage operational variables
- `STO_IN[Nodes_Storages, Hours]` - *Storage inflow* - ``STO^{\text{IN}}_{(n,sto),h}``
    - Detail: Storage energy inflow
    - Bound: lower_bound=0
    - Units: MWh per time-interval; 
- `STO_OUT[Nodes_Storages, Hours]` - *Storage outflow* - ``STO^{\text{OUT}}_{(n,sto),h}``
    - Detail: Storage energy outflow
    - Bound: lower_bound=0
    - Units: MWh per time-interval; 
- `STO_L[Nodes_Storages, Hours]` - *Storage level* - ``STO^{\text{L}}_{(n,sto),h}``
    - Detail: Storage energy level
    - Bound: lower_bound=0
    - Units: MWh at a given time-interval; 

### Generation and storage capacity variables
- `N_TECH[Nodes_Techs]` - *Technology capacity* - ``N^{\text{TECH}}_{(n,t)}``
    - Detail: Generation technology capacity 
    - Bound: lower_bound=0
    - Units: MW; 
- `N_STO_E[Nodes_Storages]` - *Storage build energy* - ``N^{\text{E}}_{(n,sto)}``
    - Detail: Storage energy capacity
    - Bound: lower_bound=0
    - Units: MWh; 
- `N_STO_P[Nodes_Storages]` - *Storage capacity* - ``N^{\text{P}}_{(n,sto)}``
    - Detail: Storage (charge/discharge) power capacity
    - Bound: lower_bound=0
    - Units: MW; 
- `N_SYNC[Nodes_SynCons]` - *SynCon capacity* - ``N^{\text{SYNC}}_{(n,s)}``
    - Detail: Synchronous condenser capacity
    - Bound: lower_bound=0
    - Units: MW; 

#### Transmission network variables
- `FLOW[Arcs,Hours]` - *Internodal flow* - ``F_{(n_F,n_T),h}``
    - Detail: Power flow between nodes in topology
    - Bound: free
    - Units: MWh; 
- `N_IC_EXP[Arcs]` - *Internodal flow expansion* - ``N^{\text{IC\_EXP}}_{(n_F,n_T)}``
    - Detail: Power transmission (flow) capacity expansion betweeen zones
    - Bound: lower_bound=0
    - Units: MW; 
- `N_REZ_EXP[REZones]` - *Renewable capacity expand* - ``N^{\text{RX}}_z``
    - Detail: Renewable technology transmission capacity built
    - Bound: lower_bound=0
    - Units: MW; 
- `N_REZ_EXP_TX[REZones]` - *REZ expand by transmission*
    - Detail: REZ transmission capacity built when transmission network is expanded
    - Bound: lower_bound=0
    - Units: MW; 

#### Hydrogen production operational variables
- `H2_P2G[Nodes_P2G, Hours]` - *H2 power to gas* - ``H2^{\text{P2G}}_{(n,p),h}``
    - Detail: Power-to-gas hydrogen production
    - Bound: lower_bound=0
    - Units: MWh per time-interval; 
- `H2_G2P[Nodes_G2P, Hours]` - *H2 gas to power* - ``H2^{\text{G2P}}_{(n,g),h}``
    - Detail: Hydrogen Gas-to-power energy conversion
    - Bound: lower_bound=0
    - Units: MWh per time-interval; 
- `H2_GS_L[Nodes_GasStorages, Hours]` - *H2 storage level* - ``H2^{\text{GS\_L}}_{(n,p),h}``
    - Detail: Current gas storage level
    - Bound: lower_bound=0
    - Units: tonne-H2 at a given time-interval; 
- `H2_GS_IN[Nodes_GasStorages, Hours]` - *H2 storage inflow*
    - Detail: Current gas storage input
    - Bound: lower_bound=0
    - Units: tonne-H2 at a given time-interval; 
- `H2_GS_OUT[Nodes_GasStorages, Hours]` - *H2 storage outflow*
    - Detail: Current gas storage output
    - Bound: lower_bound=0
    - Units: tonne-H2 at a given time-interval; 

#### Hydrogen capacity variables
- `N_P2G[Nodes_P2G]` - *H2 P2G capacity* - ``N^{\text{P2G}}_{(n,p)}``
    - Detail:  Power-to-gas capacity
    - Bound: lower_bound=0
    - Units: MW; 
- `N_G2P[Nodes_G2P]` - *H2 G2P capacity* - ``N^{\text{G2P}}_{(n,g)}``
    - Detail: Gas-to-power capacity
    - Bound: lower_bound=0
    - Units: MW; 
- `N_GS[Nodes_GasStorages]` - *H2 storage capacity* - ``N^{\text{GS}}_{(n,gs)}``
    - Detail: Gas storage capacity
    - Bound: lower_bound=0
    - Units; tonne-H2 ; 

### Objective function 

The default objective function represents the total system cost. It has an operational component and a capacity investment component. We will write
```math
Z = Z^{\text{cap}} + \sum_{h \in \text{Hours}} Z_h^{\text{op}}
```
where each cost component comprises further subcomponents, as outlined next.

#### Operational costs
For each time step index ``h \in \text{Hours}``, operational costs $Z_h^{\text{op}}$ are the sum of applicable variable costs (typically on a cost-per-MWh basis) for
-  generation, due to 
    - operation and maintenance (O&M) : ``vc``
    - fuel : ``fc``
    - carbon emissions : ``co2``
- storage charge and discharge
    - operation and maintenance (O&M) : ``vc``
- hydrogen production and hydrogen-fuelled power generation
    - operation and maintenance (O&M) : ``vc``

Mathematically,
```math
\begin{align*}
Z_h^{\text{op}} &= 
 \sum_{n \in \text{Nodes}} \sum_{t \in \text{Techs}} (vc_{(n,t)} + \eta_{(n,t)} fc_{(n,t)} + co2_{(n,t)}) G_{(n,t),h}\\
&  +  \sum_{n \in \text{Nodes}} \sum_{sto \in \text{Storages}} vc_{(n,sto)} (STO^{\text{IN}}_{(n,sto),h} + STO^{\text{OUT}}_{(n,sto),h}) \\
&  +  \sum_{n \in \text{Nodes}} \left[ \sum_{p \in \text{P2G}} vc_{(n,p)} H2^{\text{P2G}}_{(n,p),h} + \sum_{g \in \text{G2P}} vc_{(n,g)} H2^{\text{G2P}}_{(n,g),h} \right]               
\end{align*}
```
where the parameter ``\eta`` is a conversion efficiency of fuel to electrical power (storage round-trip efficiencies are handled in constraints).

In code, the objective function is
```julia
@objective(m, Min, Z)
@constraint(m, ObjectiveFunction, Z ==
sum(MarginalCost[n,t] * G[(n,t),h] for (n,t) in Nodes_Dispatch, h in Hours)

+ sum(MarginalCost[n,sto] * (STO_OUT[(n,sto),h] + STO_IN[(n,sto),h]) for (n,sto) in Nodes_Storages, h in Hours)

+ sum(MarginalCost[n,p2g] * H2_P2G[(n,p2g),h] for (n,p2g) in Nodes_P2G, h in Hours)
+ sum(MarginalCost[n,g2p] * H2_G2P[(n,g2p),h] for (n,g2p) in Nodes_G2P, h in Hours)
...
```
where `MarginalCost` encapsulates the cost components above, and the sum over all time steps is represented. Note the use of particular mapping subsets (e.g. `Nodes_Dispatch`) to refine the summation to only the terms with applicable data.

Additional operational terms for load following costs are also included as
```julia
+ LoadIncreaseCost[n,t] * sum(G_UP[(n,t),h] for (n,t) in Nodes_Dispatch, h in Hours2)
+ LoadDecreaseCost[n,t] * sum(G_DO[(n,t),h] for (n,t) in Nodes_Dispatch, h in Hours)
...
```
that serve as soft-penalty proxies for the relative flexibility of different dispatchable generation technologies.

#### Capacity investment costs 
For each modelling year, the fixed cost (annual) and the investment cost (amortised relative to the lifetime of the technology) are combined to form ``Z^{\text{cap}}``, where
```math
Z^{\text{cap}} = Z^{\text{fix}} + Z^{\text{inv}}.
```
Here
```math
\begin{align*}
Z^{\text{inv}} &= 
     \sum_{n \in \text{Nodes}} \sum_{t \in \text{Techs}} ic_{(n,t)} N^{\text{TECH}}_{(n,t)} \\
  &+ \sum_{n \in \text{Nodes}} \sum_{sto \in \text{Techs}} ic^P_{(n,sto)} N^{\text{P}}_{(n,sto)} + ic^E_{(n,sto)} N^{\text{E}}_{(n,sto)} \\
  &+ \sum_{n \in \text{Nodes}} \sum_{s \in \text{SynCons}} ic^{SYNC}_{(n,s)} N^{\text{SYNC}}_{(n,s)} \\
  &+ \sum_{n \in \text{Nodes}} \left[ \sum_{p \in \text{P2G}} ic^{P2G}_{(n,p)} N^{\text{P2G}}_{(n,p)} +  \sum_{g \in \text{G2P}} ic^{G2P}_{(n,g)} N^{\text{G2P}}_{(n,g)}
  +  \sum_{gs \in \text{GasSto}} ic^{GS}_{(n,gs)} N^{\text{GS}}_{(n,gs)}  \right] \\
  &+ \sum_{(n_F,n_T) \in \text{Arcs}} ic^{IC\_EXP}_{(n_F,n_T)} N^{\text{IC\_EXP}}_{(n_F,n_T)}
\end{align*}
```
where the parameters ``ic`` refer to the present year investment costs. The objective term ``Z^{\text{fix}}`` is analogously defined after replacing investment costs with fixed costs.


Continuing from above, the following terms are appear in the objective function as
```julia
+ sum(InvestmentCost[n,t] * N_TECH[(n,t)] for (n,t) in Nodes_Techs_New)

+ sum(InvestmentCostPower[n,sto] * N_STO_P[(n,sto)] for (n,sto) in Nodes_Storages_New)
+ sum(InvestmentCostEnergy[n,sto] * N_STO_E[(n,sto)] for (n,sto) in Nodes_Storages_New)

+ sum(InvestmentCostSynCon[n,syn] * N_SYNC[(n,syn)] for (n,syn) in Nodes_SynCons)

+ sum(InvestmentCost[n,p2g] * N_P2G[(n,p2g)] for (n,p2g) in Nodes_P2G)
+ sum(InvestmentCost[n,g2p] * N_G2P[(n,g2p)] for (n,g2p) in Nodes_G2P)
+ sum(InvestmentCost[n,gs] * N_GS[(n,gs)] for (n,gs) in Nodes_GasStorages)

+ sum(ArcWeight[from,to]*InvestmentCostTransExp[from,to] * N_IC_EXP[(from,to)] for (from,to) in Arcs)

+ sum(FixedCost[n,t] * N_TECH[(n,t)] for (n,t) in Nodes_Techs)
+ sum(FixedCost[n,sto] * N_STO_P[(n,sto)] for (n,sto) in Nodes_Storages)

+ sum(FixedCost[n,p2g] * N_P2G[(n,p2g)] for (n,p2g) in Nodes_P2G)
+ sum(FixedCost[n,g2p] * N_G2P[(n,g2p)] for (n,g2p) in Nodes_G2P)
+ sum(FixedCost[n,gs] * N_GS[(n,gs)] for (n,gs) in Nodes_GasStorages)
```

Note that while this formulation does not explicitly display discounting factors, these are applied in a separate pre-processing step.

### Constraints

The model constraints are given here in both mathematical form and code implementation form.

For compactness we abbreviate `Nodes` to ``N``, `Techs` to ``T``, `Storages` to ``S``, and `Arcs` to ``A``. 

In the code formulation, constants of the form `CapAdd[:CapacitySymbol][i,j,...]` denote capacity quantities from the optimal capacity values from the model's solution in the previous time-step. Their use is to distinguish new capacity decision variables from prior capacity expansion data. Note that in the mathematical formulation there is no separate term corresponding to `CapAdd`. Instead, the capacity variables (e.g. ``N^{\text{TECH}}``) should be read as being inclusive of prior capacity expansion values, unless otherwise stated.

#### _Energy balance and load_

*Definition of REZone generation book-keeping variables*: for each ``z \in `` `REZones`,

```math
G^{\text{REZ}}_{z,h} = 
    \sum_{t | (z,t) \in N \times T} G_{(z,t),h}
    + \sum_{sto | (z,sto) \in N \times S} (STO^{\text{IN}}_{(z,sto),h} - STO^{\text{OUT}}_{(z,sto),h}) 
```

```julia
@constraint(m, REZoneGen[rez=REZones,h=Hours],
    G_REZ[rez,h] ==
        sum(G[(z,t),h] for (z,t) in Nodes_Techs if z == rez)
        +  sum(STO_OUT[(q,sto),h] - STO_IN[(q,sto),h] for (q,sto) in Nodes_Storages if q == rez)
);
```
*Definition of TxZone generation book-keeping variables*: for each ``z \in`` `TxZones` and ``h \in`` `Hours`,

```math
G^{\text{TxZ}}_{(z,t),h} = 
    \sum_{t | (z,t) \in N \times T} G_{(n,t),h}
    + \sum_{sto | (z,sto) \in N \times S} (STO^{\text{IN}}_{(z,sto),h} - STO^{\text{OUT}}_{(z,sto),h}) \\
    + \sum_{(rez,t) \in N \times T \, | \; rez \, \uparrow z } G^{\text{REZ}}_{rez,h} - \sum_{(n_F,n_T) \in A \, | \, n_F = z } F_{(n_F,n_T),h} 
```
where ``rez \, \uparrow z``  means the renewable enegy zone ``rez`` is connected to the transmission zone ``z``. 

```julia
@constraint(m, TxZoneGen[zone=TxZones,h=Hours],
    G_TxZ[zone,h] ==
        sum(G[(z,t),h] for (z,t) in Nodes_Techs if z == zone)
        +  sum(STO_OUT[(z,sto),h] - STO_IN[(z,sto),h] for (z,sto) in Nodes_Storages if z == zone)
        +  sum(G_REZ[rez,h] for (rez, z) in Nodes_Promotes if z == zone)
        -  sum(FLOW[(from,to),h] for (from,to) in Arcs if from == zone)
);
```

*Energy balance at each demand node: supply equals or exceeds demand*: for each ``n \in `` `DemandZones` and ``h \in`` `Hours`,

```math
G^{\text{TxZ}}_{(n,t),h} +  \sum_{g \in \text{G2P}}  H2^{\text{G2P}}_{(n,g),h}
\geq 
\alpha^{\text{Tx}} * \left[ D_{n,h} +  \sum_{p \in \text{P2G}}  H2^{\text{P2G}}_{(n,p),h} \right]
```
where ``D`` is the array of hourly demands at nodes, and ``\alpha^{\text{Tx}}`` is a constant factor accounting for transmission line losses.

```julia
@constraint(m, EnergyBalance[n=DemandZones,h=Hours],
    sum(G_TxZ[zone,h] for (zone, d) in Nodes_Demand if d == n)
    + sum(H2_G2P[(zone,g2p),h] for (zone,g2p) in Nodes_G2P if zone == n)
    >=
    loss_factor_tx*(
        Load[n,h]
        + sum(H2_P2G[(zone,p2g),h] for (zone,p2g) in Nodes_P2G if zone == n))
);
```

*Energy flow reflexive constraint*: for each ``(n_F,n_T) \in`` `Arcs` and ``h \in`` `Hours`,
```math
    F_{(n_F,n_T),h} = - F_{(n_T,n_F),h}
```

```julia
@constraint(m, FlowEnergyReflex[(from,to)=Arcs,h=Hours; from in Arcs_From],
    FLOW[(from,to),h] + FLOW[(to,from),h] == 0
);
```

*Energy flow bounds*: for each ``(n_F,n_T) \in`` `Arcs` and ``h \in`` `Hours`, 
```math
    F_{(n_F,n_T),h} \leq C^{\text{Tx}}_{(n_F,n_T)} + N^{\text{IC\_EXP}}_{(n_F,n_T)}
```
```julia
@constraint(m, FlowEnergyUpperBound[(from,to)=Arcs,h=Hours],
    FLOW[(from,to),h] <= time_ratio * ( TransferCapacity[(from,to)] + N_IC_EXP[(from,to)] + CapAdd[:N_IC_EXP][(from,to)])
);
```

*Energy flow expansion symmetry*: for each ``(n_F,n_T) \in`` `Arcs`,
```math
N^{\text{IC\_EXP}}_{(n_F,n_T)} = N^{\text{IC\_EXP}}_{(n_T,n_F)}
```
```julia
@constraint(m, FlowEnergySymmetry[(from,to)=Arcs; from in Arcs_From],
    N_IC_EXP[(from,to)] == N_IC_EXP[(to,from)]
);
```

*Generation level start*: for each ``(n,t) \in `` `Nodes_Dispatch`,
```math
    G_{(n,t),1} = G^{\uparrow}_{(n,t),1}
```
```julia
@constraint(m, GenLevelStart[(n,t)=Nodes_Dispatch],
            G[(n,t),Hours[1]] == G_UP[(n,t),Hours[1]]
);
```

*Generation level dynamics*: for each ``(n,t) \in `` `Nodes_Dispatch`, and for ``h \in`` `Hours`, ``h \neq 1``,
```math
    G_{(n,t),h} = G_{(n,t),h-1} + G^{\uparrow}_{(n,t),h} - G^{\downarrow}_{(n,t),h}
```
```julia
@constraint(m, GenLevelUpdate[(n,t)=Nodes_Dispatch,h=Hours2],
    G[(n,t),h] == G[(n,t),h-1] + G_UP[(n,t),h] - G_DO[(n,t),h]
);
```

*Variable upper bound on dispatchable generation by capacity*: for each ``(n,t) \in `` `Nodes_Dispatch`, and for ``h \in`` `Hours`,
```math
    G_{(n,t),h} \leq C^{\text{Derate}}_{(n,t),h} N^{\text{TECH}}_{(n,t)}
```
where ``C^{\text{Derate}}`` denotes capacity derating depending on generation technology type and hour characteristics (such as seasonal temperature).
```julia
@constraint(m, MaxGenerationDisp[(n,t)=Nodes_Dispatch,h=Hours],
    G[(n,t),h] <= CapacityDerating[n,t,h] * time_ratio * (N_TECH[(n,t)] + CapAdd[:N_TECH][(n,t)])
);
```

*Variable upper bound on non-dispatchable generation by capacity*: for each ``(n,t) \in `` `Nodes_Avail_Techs`, and for ``h \in`` `Hours`,
```math
G_{(n,t),h} \leq C^{\text{Avail}}_{(n,t),h} N^{\text{TECH}}_{(n,t)}
```
where ``C^{\text{Avail}}`` denotes the availability of a variable renewable source (a value between 0 and 1 inclusive).
```julia
@constraint(m, MaxGenerationNonDisp[(n,t)=Nodes_Avail_Techs,h=Hours],
    G[(n,t),h] <= Availability[n,t,h] * time_ratio * (N_TECH[(n,t)] + CapAdd[:N_TECH][(n,t)])
);
```

*Maximum capacity allowed*: for each ``(n,t) \in `` `Nodes_Techs`,
```math
    N^{\text{TECH}}_{(n,t)} \leq \widehat{N}_{(n,t)}
```
```julia
@constraint(m, MaxCapacityBound[(n,t)=Nodes_Techs; !(MaxCapacity[n,t] |> ismissing)],
    N_TECH[(n,t)] + CapAdd[:N_TECH][(n,t)] <= MaxCapacity[n,t]
);
```

*Maximum generated energy allowed*: for each ``(n,t) \in `` `Nodes_Techs`,
```math
    \sum_{h \in Hours} G_{(n,t),h} \leq \widehat{E}_{(n,t)}
```
```julia
@constraint(m, MaxEnergyGenerated[(n,t)=Nodes_Techs; !(MaxEnergy[n,t] |> ismissing)],
    sum(G[(n,t),h] for h in Hours) <= MaxEnergy[n,t]
);
```

*Renewable energy zone build limits*: for each ``z \in `` `REZones` and ``h \in`` `Hours`,
```math
\sum_{h \in Hours} G_{(z,t),h} \leq \bar{C}^{\text{RX}}_z + N^{\text{RX}}_z
```
where the ``\bar{C}^{\text{RX}}`` term denotes the available capacity to build in a zone before further expansion capacity ``N^{\text{RX}}`` is required.
```julia
@constraint(m, REZBuildLimits[rez=REZones,h=Hours],
    sum(G[(z,t),h] for (z,t) in Nodes_Avail_Techs if z == rez)
        <= time_ratio*(TotalBuildCap[rez] +  N_REZ_EXP[rez] + CapAdd[:N_REZ_EXP][rez])
);
```

*Renewable energy zone expansion limits*: for each ``z \in `` `REZones`,
```math
N^{\text{RX}}_z \leq \widehat{N}^{\text{RX}}_z 
```
```julia
@constraint(m, REZExpansionBound[rez=REZones],
    N_REZ_EXP[rez] + trunc(CapAdd[:N_REZ_EXP][rez]) <= ExpansionLimit_REZ[rez]
);
```

<!-- *Renewable energy zone expansion link to transmission expansion*: -->
