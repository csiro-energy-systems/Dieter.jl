-- cf. ../base/technologies.sql for analogue for general non-storage technology
SELECT
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
OvernightCost as H2OvernightCost,
-- OvernightCostEnergy,
CurtailmentCost,
LoadIncreaseCost,
LoadDecreaseCost,
InvestLifetime AS Lifetime,
--InvestRecovery,
MaxCapacity,
-- MaxEnergy,
StartLevel,
-- MaxEnergyToPowerRatio,
IsRenewable AS Renewable,
IsDispatchable AS Dispatchable,
IsStorage AS Storage,
IsHybrid AS Hybrid
FROM
H2_Tech
-- WHERE FuelType = 'Hydrogen'
