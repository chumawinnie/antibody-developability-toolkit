# Auto-generated PyMOL script
# Visualizes biophysical liabilities on trastuzumab Fab

fetch 1n8z
remove resn HOH
bg_color white
color grey80, chain A or chain B
color palegreen, chain C

select aggregation, (chain A and resi 8) or (chain A and resi 11) or (chain A and resi 15) or (chain A and resi 40) or (chain A and resi 53) or (chain A and resi 54) or (chain A and resi 58) or (chain A and resi 59) or (chain A and resi 80) or (chain A and resi 95) or (chain A and resi 98) or (chain A and resi 110) or (chain A and resi 112) or (chain A and resi 119) or (chain A and resi 125) or (chain A and resi 141) or (chain A and resi 146) or (chain A and resi 153) or (chain A and resi 154) or (chain A and resi 163) or (chain A and resi 181) or (chain A and resi 184) or (chain A and resi 191) or (chain A and resi 201) or (chain A and resi 204) or (chain A and resi 205) or (chain A and resi 214) or (chain B and resi 2) or (chain B and resi 5) or (chain B and resi 11) or (chain B and resi 14) or (chain B and resi 18) or (chain B and resi 23) or (chain B and resi 27) or (chain B and resi 40) or (chain B and resi 41) or (chain B and resi 72) or (chain B and resi 88) or (chain B and resi 93) or (chain B and resi 104) or (chain B and resi 110) or (chain B and resi 115) or (chain B and resi 121) or (chain B and resi 130) or (chain B and resi 132) or (chain B and resi 153) or (chain B and resi 156) or (chain B and resi 159) or (chain B and resi 165) or (chain B and resi 166) or (chain B and resi 170) or (chain B and resi 174) or (chain B and resi 175) or (chain B and resi 176) or (chain B and resi 177) or (chain B and resi 182) or (chain B and resi 192) or (chain B and resi 196) or (chain B and resi 202) or (chain B and resi 209) or (chain B and resi 214) or (chain B and resi 220)
color orange, aggregation
show sticks, aggregation and not name N+C+O

select deamidation, (chain A and resi 158) or (chain B and resi 55) or (chain B and resi 84) or (chain B and resi 162)
color red, deamidation
show sticks, deamidation

select oxidation, (chain B and resi 110)
color magenta, oxidation
show sticks, oxidation

zoom chain A or chain B
ray 1600, 1200
png trastuzumab_liabilities.png
save trastuzumab_liabilities_session.pse
