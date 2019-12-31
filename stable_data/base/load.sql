-- SELECT
-- Demand_Series.TimeIndex,
-- Demand_Series.RegionID,
-- Demand_Series.Demand AS Load
-- FROM Demand_Series
SELECT
Demand_Series.TimeIndex,
replace(RegionID,'DE','Demand_DE') AS DemandRegion,
-- replace(replace(RegionID,'DE','Demand_DE'),'FR','Demand_FR') AS DemandRegion,
Demand_Series.Demand AS Load
FROM Demand_Series
UNION
SELECT
Demand_Series.TimeIndex,
replace(RegionID,'DE','Demand_FR') AS DemandRegion, 
round(Demand_Series.Demand*64.4/80.7,1) AS Load   -- 64.4M:80.7M FR:DE pop. (est. 2015)
FROM Demand_Series