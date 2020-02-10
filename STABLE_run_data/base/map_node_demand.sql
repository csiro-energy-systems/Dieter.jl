SELECT
RegionID as Nodes,
DemandRegion,
IncludeLevel as IncludeFlag
FROM Regions
WHERE DemandRegion IS NOT NULL
AND IncludeLevel=1
