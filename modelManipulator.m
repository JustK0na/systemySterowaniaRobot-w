clearvars;
syms Q11 Q12 Q13 Q22 Q23 Q33;
syms l1 l2 l3 m1 m2 m3
syms q1 q2 q3
syms dq1 dq2 dq3
syms ddq1 ddq2 ddq3
syms g

lenghts = [0.5 0.5 0.5];
masses = [10 2 1];

Q11 = 1/3 * m1 * l1^2 + m2*l1^2 + 1/3 * m2 * l2^2 + m2*l1*l2*cos(q2) + 1/3 * m3 *l3^2 + m3*(l1^2 + l2^2 + l2*l3*cos(q3) + 2*l1*l2*cos(q2)+l1*l3*cos(q2+q3))  ;
Q12 = 1/3 * m2 * l2^2 + 1/2 * m2 * l1 * l2 * cos(q2) + 1/3 * m3 * l3^2 + m3*(l2^2+l2*l3*cos(q3)+l1*l2*cos(q2)+1/2*l1*l3*cos(q2+q3));
Q13 = 1/3 * m3 * l3^2 + 1/2 * m3 * (l2*l3*cos(q3) + l1*l3*cos(q2+q3));
Q22 = 1/3 * m2 * l2^2 + 1/3 * m3 * l3^2 + m3*(l2*l3*cos(q3) + l2^2);
Q23 = 1/3 * m3 * l3^2 + 1/2 * m3 * l2 * l3 * cos(q3);
Q33 = 1/3 * m3 * l3^2;
Msym = [Q11 Q12 Q13;
     Q12 Q22 Q23;
     Q13 Q23 Q33];

C11 = -1/2 * dq2 * (m2*l1*l2*sin(q2) + 2*m3*l1*l2*sin(q2) + m3*l1*l3*sin(q2+q3)) - 1/2*dq3*(m3*l2*l3*sin(q3) + m3*l1*l3*sin(q2+q3));
C12 = -1/2 * (dq1 + dq2) * (m2*l1*l2*sin(q2) + 2*m3*l1*l2*sin(q2) + m3*l1*l3*sin(q2+q3)) - 1/2 * dq3 * (m3*l2*l3*sin(q3) + m3*l1*l3*sin(q2+q3));
C13 = -1/2 * (dq1 + dq2 + dq3) * (m3*l2*l3*sin(q3) + m3*l1*l3*sin(q2+q3));
C21 = 1/2 * dq1 * (m2*l1*l2*sin(q2) + 2*m3*l1*l2*sin(q2) + m3*l1*l3*sin(q2+q3)) - 1/2*dq3*m3*l2*l3*sin(q3);
C22 = -1/2 * dq3 * m3 * l2 * l3 * sin(q3);
C23 = -1/2 *(dq1 + dq2 + dq3)*m3*l2*l3*sin(q3);
C31 = 1/2 * dq1 * (m3*l2*l3*sin(q3) + m3*l1*l3*sin(q2+q3)) + 1/2 * dq2*m3*l2*l3*sin(q3);
C32 = 1/2 * (dq1 + dq2)*m3*l2*l3*sin(q3);

Csym = [C11 C12 C13;
    C21 C22 C23;
    C31 C32 0];

D1 = 1/2 * 9.81 * ( (m1 + 2*m2 + 2*m3)*l1*cos(q1) + (m2 + 2*m3)*l2*cos(q1+q2) + m3*l3*cos(q1+q2+q3));
D2 = 1/2 * 9.81 * ( (m2 + 2*m3)*l2*cos(q1+q2) + m3*l3*cos(q1+q2+q3));
D3 = 1/2 * 9.81 * m3*l3*cos(q1+q2+q3);

Dsym = [D1;
    D2;
    D3]; 


steps = 0:0.1:(2*pi)
numSteps = length(steps)

Eigen = zeros(3, numSteps^2);
idx = 1;
qs = [0, 0, 0]

M = subs(Msym, [l1 l2 l3 m1 m2 m3], [lenghts, masses])
C = subs(Csym, [l1 l2 l3 m1 m2 m3], [lenghts, masses])
D = double(subs(Dsym, [l1 l2 l3 m1 m2 m3 q1 q2 q3], [lenghts, masses, qs]))
% 
% %=0:0.1:(2*pi)
%     for j=0:0.1:(2*pi)
%         for k=0:0.1:(2*pi)
%             qs = [j k];
%             M = subs(Msym, [l1 l2 l3 m1 m2 m3 q2 q3], [lenghts, masses, qs]);
%             Eigen(:,idx) = eigs(double(M)); 
%             idx=idx+1;
%           %  fprintf('j: %d, k: %d \n', j, k);
%         end
%     end
% Zdefiniowanie wektorów q i dq
q = [q1; q2; q3];
dq = [dq1; dq2; dq3];

% Obliczenie pochodnej macierzy bezwładności M po czasie (dM)
dM = zeros(3,3, 'sym');
for i = 1:3
    dM = dM + diff(Msym, q(i)) * dq(i);
end

test_symetrii = simplify(dM - (Csym + Csym.'));

disp('Wynik testu dM = C + C^T (powinna być macierz zerowa):');
disp(test_symetrii);