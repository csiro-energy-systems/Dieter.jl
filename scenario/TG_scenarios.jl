
# %% Mappings from Scenario designations to parameters
Scenario_Number_Dict = Dict(
      "Scen1_BAU" => 1,  # Current Trends
      "Scen2_DDC" => 2,  # Deep Decarbonisation
      "Scen3_PRP" => 3,  # Prosumer Power
      "Scen4_DID" => 4,  # De-Industrialisation Death Spiral
      "Scen5_BTF" => 5,  # Back to the Future
      "Scen6_NGA" => 6,  # NSW Goes It Alone
      "Scen7_DCU" => 7,  # The DC Universe
      "Scen8_LCE" => 8   # Australia's Low Cost Energy Advantage
)

Scenario_TxCost_Scaling_Dict = Dict(
      "Scen1_BAU" => 1,  # Current Trends
      "Scen2_DDC" => 1,  # Deep Decarbonisation
      "Scen3_PRP" => 1,  # Prosumer Power
      "Scen4_DID" => 1,  # De-Industrialisation Death Spiral
      "Scen5_BTF" => 1,  # Back to the Future
      "Scen6_NGA" => 1,  # NSW Goes It Alone
      "Scen7_DCU" => 1,  # The DC Universe
      "Scen8_LCE" => 1   # Australia's Low Cost Energy Advantage
)

Scen_ISP_Map = Dict("Scen1_BAU" => "Central", "Scen2_DDC" => "HighVRE")

Scen_FuelCost_Map = Dict("Scen1_BAU" => "Neutral", "Scen2_DDC" => "Fast")

Scenario_BattVPP_Dict = Dict(
      "Scen1_BAU" => 2,  # Current Trends
      "Scen2_DDC" => 3,  # Deep Decarbonisation
      "Scen3_PRP" => 4,  # Prosumer Power
      "Scen4_DID" => 1,  # De-Industrialisation Death Spiral
      "Scen5_BTF" => 1,  # Back to the Future
      "Scen6_NGA" => 3,  # NSW Goes It Alone
      "Scen7_DCU" => 2,  # The DC Universe
      "Scen8_LCE" => 3   # Australia's Low Cost Energy Advantage
)

# A value of `missing` means H2 tech. not included in model,
# 0 means H2 tech. is included
#    and there is a simple minimum amount of hydrogen in a year from Electrolysers.
# 2 means H2 tech. included
#    and there is a round-trip storage path of Electrolysers, GasStorage and Gas-To-Power.
Scenario_H2_Map = Dict(
      "Scen1_BAU" => 0,  # Current Trends
      "Scen2_DDC" => 0,  # Deep Decarbonisation
      "Scen3_PRP" => 0,  # Prosumer Power
      "Scen4_DID" => 0,  # De-Industrialisation Death Spiral
      "Scen5_BTF" => 0,  # Back to the Future
      "Scen6_NGA" => 0,  # NSW Goes It Alone
      "Scen7_DCU" => 0,  # The DC Universe
      "Scen8_LCE" => 0   # Australia's Low Cost Energy Advantage
)

# %% Allowed sets for scenarios

Allowed_Scenarios = keys(Scenario_Number_Dict)
Allowed_Years = [2020, 2040, 2050]

# %% Technology type definitions
scen_types = Dict{Symbol, String}()
BattEnergyType = "N_BattEnergy"
HydPumpEnergyType = "N_HydPumpEnergy"

H2ElectrolyserType = "N_Electrolyser"
H2RecipEngType = "N_RecipH2"

scen_types[:BattEnergyType] = BattEnergyType
scen_types[:HydPumpEnergyType] = HydPumpEnergyType
scen_types[:H2ElectrolyserType] = H2ElectrolyserType
scen_types[:H2RecipEngType] = H2RecipEngType
