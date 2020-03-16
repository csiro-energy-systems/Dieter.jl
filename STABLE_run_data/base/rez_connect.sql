SELECT
REZ_Build.RezID AS Region,
TechID AS Technologies,
ConnectCostMW AS ConnectCost
FROM
REZ_Build
JOIN
(
SELECT
TechID, RegionID
FROM
Node_Tech_All_Map
WHERE
Status = 'NewEntrant'
AND IncludeFlag = 1
AND TechTypeCategory = 'Renewable'
)
ON (REZ_Build.RezID = RegionID)
