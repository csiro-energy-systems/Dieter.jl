SELECT
RegionID as Nodes,
DemandRegion,
1 as IncludeFlag
FROM Regions
WHERE DemandRegion IS NOT NULL