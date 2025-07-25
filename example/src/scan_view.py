import chimera
from chimera import runCommand

runCommand('sel #0.2- & protein')
runCommand('~ribbon sel')
runCommand('~disp sel')
runCommand('sel #0.1 & protein')
runCommand('color pink,r sel')
runCommand('sel #0.1 & protein & element.C')
runCommand('color pink,a sel')
runCommand('sel #0.2-9999 & :FE1')
runCommand('~disp sel')
runCommand('sel #0.1')
runCommand('transp 96,r sel')
runCommand('sel #0.1-69 & :HPD & element.C')
runCommand('color #bfffbf sel')
runCommand('sel #0.70 & :HPD & element.C')
runCommand('color light blue sel')
runCommand('sel #0.2-69 & :HPD')
runCommand('transp 75,a sel')
runCommand('sel #0.2-69 & :OH1')
runCommand('transp 75,a sel')
runCommand('sel #0.1,70:HPD@h14,h15,h18,h19,h23,h31')
runCommand('disp sel')
runCommand('background solid white')
runCommand('sel #0.1 & protein')
runCommand('~disp sel')
runCommand('sel solvent')
runCommand('~disp sel')
runCommand('save <PATH>')











