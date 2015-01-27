bg_color white
hide all
show cartoon, chain A
set cartoon_transparency, 0.7
color gray85, chain A
select heteroatoms, (hetatm and not resn HOH) AND chain A
show sticks, heteroatoms
color magenta, heteroatoms
select LYS2, resi 2 AND chain A
select aroundLYS2, resi 3+22+53+54+55+499+ AND chain A
select LYS72, resi 72 AND chain A
select aroundLYS72, resi 71+73+398+ AND chain A
select LYS95, resi 95 AND chain A
select aroundLYS95, resi 26+27+28+91+92+94+96+345+346+349+427+436+445+481+482+494+496+ AND chain A
select LYS99, resi 99 AND chain A
select aroundLYS99, resi 94+97+98+100+103+104+ AND chain A
select LYS120, resi 120 AND chain A
select aroundLYS120, resi 118+121+122+123+131+315+318+319+403+476+ AND chain A
select LYS133, resi 133 AND chain A
select aroundLYS133, resi 129+130+134+136+173+ AND chain A
select LYS145, resi 145 AND chain A
select aroundLYS145, resi 115+116+117+144+146+147+181+182+183+207+209+256+311+371+437+484+ AND chain A
select LYS198, resi 198 AND chain A
select aroundLYS198, resi 194+195+197+199+201+226+227+483+511+ AND chain A
select LYS221, resi 221 AND chain A
select aroundLYS221, resi 217+218+220+222+225+250+251+ AND chain A
select LYS308, resi 308 AND chain A
select aroundLYS308, resi 113+206+230+307+309+401+405+409+419+446+447+448+459+ AND chain A
select LYS362, resi 362 AND chain A
select aroundLYS362, resi 8+41+44+361+363+384+417+443+507+ AND chain A
show sticks, LYS* AND chain A
color red, LYS* AND chain A
color tv_blue, around* AND chain A
select histidines, resn his AND chain A
show sticks, histidines
show spheres, heteroatoms
deselect
