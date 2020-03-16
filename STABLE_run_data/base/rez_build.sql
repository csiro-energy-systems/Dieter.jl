SELECT
REZ_Build.RezID AS Region,
WindLimit,
SolarLimit,
TotalBuildCap,
ExpansionLimit,
TransExpansionCostMW AS TransExpansionCost
-- ConnectCostMW AS ConnectCost
FROM REZ_Build
