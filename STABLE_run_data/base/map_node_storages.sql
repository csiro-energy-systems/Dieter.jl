SELECT
RegionID as Nodes,
TechID as Technologies,
IncludeFlag
FROM
Node_Tech_All_Map
WHERE Node_Tech_All_Map.IsStorage=1
--
-- SELECT
-- RegionID as Nodes,
-- TechID as Technologies,
-- 1 as IncludeFlag
-- FROM
-- Tech_All_Parameters
-- WHERE Tech_All_Parameters.IsStorage=1
--
-- SELECT
-- Map_Region_Tech.RegionID as Nodes,
-- Map_Region_Tech.TechTypeID as Technologies,
-- Map_Region_Tech.IncludeFlag
-- FROM
-- Map_Region_Tech JOIN ...
