bg_color white
hide all
show cartoon, chain B
set cartoon_transparency, 0.7
color gray85, chain B
select heteroatoms, (hetatm and not resn HOH) AND chain B
show sticks, heteroatoms
color magenta, heteroatoms
select LYS2, resi 2 AND chain B
select aroundLYS2, resi 3+22+53+54+55+398+ AND chain B
select LYS72, resi 72 AND chain B
select aroundLYS72, resi 70+71+73+456+ AND chain B
select LYS95, resi 95 AND chain B
select aroundLYS95, resi 26+27+28+91+92+94+96+345+346+349+382+405+453+ AND chain B
select LYS99, resi 99 AND chain B
select aroundLYS99, resi 94+97+98+100+103+104+ AND chain B
select LYS120, resi 120 AND chain B
select aroundLYS120, resi 118+121+122+123+131+315+318+319+ AND chain B
select LYS133, resi 133 AND chain B
select aroundLYS133, resi 129+130+134+136+173+ AND chain B
select LYS145, resi 145 AND chain B
select aroundLYS145, resi 115+116+117+118+144+146+147+181+182+183+207+209+256+311+371+420+425+ AND chain B
select LYS198, resi 198 AND chain B
select aroundLYS198, resi 194+195+197+199+201+226+227+ AND chain B
select LYS221, resi 221 AND chain B
select aroundLYS221, resi 217+218+220+222+225+250+251+463+ AND chain B
select LYS308, resi 308 AND chain B
select aroundLYS308, resi 113+206+230+307+309+417+424+428+435+444+448+470+503+513+ AND chain B
select LYS362, resi 362 AND chain B
select aroundLYS362, resi 8+41+44+361+363+386+393+419+ AND chain B
show sticks, LYS* AND chain B
color red, LYS* AND chain B
color tv_blue, around* AND chain B
select histidines, resn his AND chain B
show sticks, histidines
show spheres, heteroatoms
deselect
