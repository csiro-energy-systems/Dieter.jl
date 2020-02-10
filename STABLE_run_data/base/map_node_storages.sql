-- SELECT
-- Map_Region_Tech.RegionID as Nodes,
-- Map_Region_Tech.TechTypeID as Technologies,
-- Map_Region_Tech.IncludeFlag
-- FROM
-- Map_Region_Tech JOIN ...
SELECT
RegionID as Nodes,
TechID as Technologies,
1 as IncludeFlag
FROM
Tech_All_Parameters
WHERE Tech_All_Parameters.IsStorage=1
