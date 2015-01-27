bg_color white
hide all
show cartoon, chain C
set cartoon_transparency, 0.7
color gray85, chain C
select heteroatoms, (hetatm and not resn HOH) AND chain C
show sticks, heteroatoms
color magenta, heteroatoms
select LYS2, resi 2 AND chain C
select aroundLYS2, resi 3+22+53+54+55+467+ AND chain C
select LYS72, resi 72 AND chain C
select aroundLYS72, resi 70+71+73+417+476+ AND chain C
select LYS95, resi 95 AND chain C
select aroundLYS95, resi 26+27+28+91+92+94+96+345+346+349+430+444+452+455+ AND chain C
select LYS99, resi 99 AND chain C
select aroundLYS99, resi 94+97+98+100+103+104+ AND chain C
select LYS120, resi 120 AND chain C
select aroundLYS120, resi 118+121+122+123+131+315+318+319+480+493+ AND chain C
select LYS133, resi 133 AND chain C
select aroundLYS133, resi 129+130+134+136+137+173+ AND chain C
select LYS145, resi 145 AND chain C
select aroundLYS145, resi 115+116+117+118+144+146+147+181+182+183+207+209+256+311+371+393+412+ AND chain C
select LYS198, resi 198 AND chain C
select aroundLYS198, resi 194+195+197+199+201+226+227+423+485+ AND chain C
select LYS221, resi 221 AND chain C
select aroundLYS221, resi 217+218+220+222+225+250+251+ AND chain C
select LYS308, resi 308 AND chain C
select aroundLYS308, resi 113+206+230+307+309+369+388+389+406+422+435+462+483+510+ AND chain C
select LYS362, resi 362 AND chain C
select aroundLYS362, resi 8+40+41+44+361+363+381+400+416+468+ AND chain C
show sticks, LYS* AND chain C
color red, LYS* AND chain C
color tv_blue, around* AND chain C
select histidines, resn his AND chain C
show sticks, histidines
show spheres, heteroatoms
deselect
