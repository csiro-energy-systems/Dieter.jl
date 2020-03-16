SELECT
TechID AS Technologies,
Regions.RegionID AS Region,
ConnectCost
FROM
TxZ_Connect_Cost JOIN Regions
ON (Regions.DemandRegion = TxZ_Connect_Cost.DemandRegion)
WHERE RegionType = 'TxZone'
