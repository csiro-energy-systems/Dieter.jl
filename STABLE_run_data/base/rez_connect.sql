SELECT
  REZ_Build.RezID AS Region,
  TechID AS Technologies,
  ConnectCost*1000 AS ConnectCost -- conversion from $/kW to $/MW
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
