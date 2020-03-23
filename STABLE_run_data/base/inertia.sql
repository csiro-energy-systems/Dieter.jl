SELECT TechID, InertialCoeffExi AS InertialCoeff
FROM
(Tech_Map JOIN Technology_Types USING (TechTypeID))
JOIN
Tech_Inertia USING (TechTypeID)
WHERE Status LIKE '%existing%'
UNION
SELECT TechID, InertialCoeffNew AS InertialCoeff
FROM
(Tech_Map JOIN Technology_Types USING (TechTypeID))
JOIN
Tech_Inertia USING (TechTypeID)
WHERE Status LIKE '%new%'
