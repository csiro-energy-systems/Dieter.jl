-- SELECT
-- Map_Region_Tech.RegionID as Nodes,
-- Map_Region_Tech.TechTypeID as Technologies,
-- Map_Region_Tech.IncludeFlag
-- FROM
-- Map_Region_Tech JOIN ...
SELECT
RegionID as Nodes,
TechTypeID as Technologies,
1 as IncludeFlag
FROM
Technologies
JOIN Technology_Types USING (TechTypeID)
WHERE Technology_Types.IsStorage=1
