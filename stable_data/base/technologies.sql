SELECT
Technologies.TechTypeID AS Technologies,
Technologies.RegionID AS Region,
Technologies.CapacityFactor,
Technologies.FuelCost,
Technologies.FixedCost,
Technologies.VariableCost,
Technologies.OvernightCostPower AS OvernightCost,
--Technologies.OvernightCostEnergy,
Technologies.CurtailmentCost,
Technologies.LoadIncreaseCost,
Technologies.LoadDecreaseCost,
Technologies.InvestLifetime AS Lifetime,
--Technologies.InvestRecovery,
Technologies.MaxCapacity AS MaxInstallable,
Technologies.MaxEnergy,
--Technologies.StartLevel,
--Technologies.MaxEnergyToPowerRatio,
--
Technology_Types.TechTypeCategory,
Technology_Types.IsRenewable AS Renewable,
Technology_Types.IsDispatchable AS Dispatchable,
-- Technology_Types.IsStorage,
Technology_Types.FuelSource,
Technology_Types.Efficiency,
Technology_Types.CarbonContent,
CO2price AS CO2_price
FROM 
(
Technologies JOIN Technology_Types USING (TechTypeID)
)
JOIN Regions USING (RegionID) -- Regions JOIN provides CO2price
WHERE Technology_Types.IsStorage=0