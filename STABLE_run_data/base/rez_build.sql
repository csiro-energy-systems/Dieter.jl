SELECT
REZ_Build.RezID AS Nodes,
WindLimit,
SolarLimit,
TotalBuildCap,
TransExpansionCostMW AS TransExpansionCost,
ConnectCostMW AS ConnectCost
FROM REZ_Build
