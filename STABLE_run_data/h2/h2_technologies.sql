-- cf. ../base/technologies.sql for analogue for general non-storage technology
SELECT
RegionID AS Region,
TechID AS H2Technologies,
Status,
TechTypeDescriptor AS TechType,
FuelType,
TechTypeCategory AS H2Type,
-- CapacityExist AS ExistingCapacity,
FuelCost,
Efficiency,
H2Conversion,
FixedCost,
VariableCost,
OvernightCost*1000 as H2OvernightCost, -- conversion from $/kW to $/MW for energy
-- OvernightCostEnergy,
CurtailmentCost,
LoadIncreaseCost,
LoadDecreaseCost,
-- InvestLifetime AS Lifetime,
InvestRecovery AS Lifetime,
MaxCapacity,
-- MaxEnergy,
StartLevel,
-- MaxEnergyToPowerRatio,
IsRenewable AS Renewable,
IsDispatchable AS Dispatchable,
IsStorage AS Storage,
IsHybrid AS Hybrid
FROM
H2_Tech JOIN Regions
WHERE RegionType = 'TxZone'
-- WHERE FuelType = 'Hydrogen'
