-- SELECT
-- Demand_Series.TimeIndex,
-- Demand_Series.RegionID,
-- Demand_Series.Demand AS Load
-- FROM Demand_Series
SELECT
Demand_Series.TimeIndex,
replace(replace(RegionID,'DE','Demand_DE'),'FR','Demand_FR') AS DemandRegion,
Demand_Series.Demand AS Load
FROM Demand_Series