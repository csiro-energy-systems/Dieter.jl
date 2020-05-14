SELECT
-- * FROM Nodes_Def
RegionID as Nodes,
RegionType as NodeType,
SuperRegion as NodePromote,
DemandRegion,
IncludeLevel as IncludeFlag
FROM Regions
WHERE RegionType NOT IN ('RootNode','DemandRegion')