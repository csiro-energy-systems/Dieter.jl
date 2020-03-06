SELECT
RegionID as Nodes,
TechID as Technologies,
IncludeFlag
FROM
Node_Tech_All_Map
WHERE Node_Tech_All_Map.IsStorage=0
--
-- SELECT
-- RegionID as Nodes,
-- TechID as Technologies,
-- 1 as IncludeFlag
-- FROM
-- Tech_All_Parameters
-- WHERE Tech_All_Parameters.IsStorage=0
--
-- SELECT
-- Map_Region_Tech.RegionID as Nodes,
-- Map_Region_Tech.TechTypeID as Technologies,
-- Map_Region_Tech.IncludeFlag
-- FROM
-- Map_Region_Tech
-- JOIN Technology_Types USING (TechTypeID)
-- WHERE Technology_Types.IsStorage=0
