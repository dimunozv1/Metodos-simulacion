set term pdf
set out 'cfigura.pdf'
set title 'Modos normales' 
set xlabel 'r'
set ylabel 'R(r,lambda)'
plot 'lambda0.dat' w l t 'R0(r,2.40527)',  'lambda1.dat' w l t 'R0(r,5.52245)','lambda2.dat' w l t 'R0(r,8.65954)', 'lambda3.dat' w l t 'R0(r,11.8022)'
