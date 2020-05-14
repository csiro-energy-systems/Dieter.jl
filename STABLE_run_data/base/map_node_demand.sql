SELECT
RegionID AS Nodes,
RegionID AS DemandZone,  -- TxZones AS Demand Regions
-- DemandRegion,
IncludeLevel AS IncludeFlag
FROM Regions
WHERE RegionType = 'TxZone'  -- TxZones AS Demand Regions
-- WHERE DemandRegion IS NOT NULL
AND IncludeLevel=1
