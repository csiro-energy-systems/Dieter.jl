-- This query creates arcs between all nodes by default:
-- CREATE VIEW Arcs_Def AS
-- SELECT FromNodes, ToNodes, 1 as IncludeFlag
-- FROM
-- (SELECT RegionID as FromNodes, RegionType as FromType
-- FROM Regions
-- WHERE RegionType=='TxZone')
-- JOIN
-- (SELECT RegionID as ToNodes, RegionType as ToType
-- FROM Regions
-- WHERE RegionType=='TxZone')
-- ON  FromType = ToType
-- WHERE FromNodes != ToNodes
--
SELECT FromNode, ToNode, TransferCapacity, TxZoneExpansionCost
FROM
(SELECT
    FromNode,
    ToNode,
    ForwardLimit AS TransferCapacity,
    TransExpansionCost*1000 AS TxZoneExpansionCost -- conversion from $/kW to $/MW
FROM Arcs
UNION
SELECT
    ToNode AS FromNode,
    FromNode AS ToNode,
    ReverseLimit AS TransferCapacity,
    TransExpansionCost*1000 AS TxZoneExpansionCost -- conversion from $/kW to $/MW
FROM Arcs)
-- SELECT * FROM Arcs_Def
