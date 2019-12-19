SELECT
<<<<<<< HEAD
Technologies.RegionID AS Region,
Technologies.TechTypeID AS Technologies,
=======
Technologies.TechTypeID AS Technologies,
Technologies.RegionID AS Region,
>>>>>>> 51849720cbbe86ef4144d8e58f86d10de5bea383
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
<<<<<<< HEAD
LEFT JOIN Region_Scalars USING (RegionID) -- Region_Scalar JOIN provides CO2price
=======
JOIN Regions USING (RegionID) -- Regions JOIN provides CO2price
>>>>>>> 51849720cbbe86ef4144d8e58f86d10de5bea383
WHERE Technology_Types.IsStorage=0