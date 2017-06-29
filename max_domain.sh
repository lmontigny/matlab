#!/bin/bash
# 1st entry = name of the file
# 2st entry = security coef
# Usage: ./max_domain.sh myFileToRead.txt 0

octave --silent << EOF
disp("------- LOADING FILE ------")
file = load('$1');
file = file.*1e-3;
dx=max(abs(file(2,:)-file(1,:)));
security_coef=$2*dx;

a=min(file);
a(:,1) = [];
a=a-security_coef;

b=max(file);
b(:,1) = [];
b=b+security_coef;

disp("");
g=sprintf('%f ', a);
fprintf('min_domain \"%s\"\n', g)
g=sprintf('%f ', b);
fprintf('max_domain \"%s\"\n', g)

quit
EOF

#octave -nodesktop -nosplash -nodisplay -r "run ./$filename ; quit;" 
