SELECT
Technologies.TechTypeID AS Storages,
Technologies.RegionID AS Region,
Technologies.CapacityFactor,
Technologies.FuelCost,
Technologies.FixedCost,
Technologies.VariableCost,
Technologies.OvernightCostPower,
Technologies.OvernightCostEnergy,
Technologies.CurtailmentCost,
--Technologies.LoadIncreaseCost,
--Technologies.LoadDecreaseCost,
-- Technologies.InvestLifetime AS Lifetime,
Technologies.InvestRecovery AS Lifetime,
Technologies.MaxCapacity,
Technologies.MaxEnergy,
Technologies.StartLevel,
Technologies.MaxEnergyToPowerRatio,
--Technology_Types.TechTypeCategory,
--Technology_Types.IsRenewable,
Technology_Types.IsDispatchable AS Dispatchable,
Technology_Types.FuelSource,
Technology_Types.Efficiency
--Technology_Types.CarbonContent
FROM
Technologies JOIN Technology_Types USING (TechTypeID)
WHERE IsStorage=1
