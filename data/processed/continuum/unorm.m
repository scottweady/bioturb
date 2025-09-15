
close all

hold on
filename = 'timeStepInfo.txt';
tc = 2.3453; % characteristic time scale

% Get contents
fileContents = readlines(filename);

% Access time and velocity
T_string = fileContents(1 : 6 : end);
U_string = fileContents(3 : 6 : end);

Nt = length(T_string);

U_norm = zeros(Nt - 1, 2);

for nt = 1 : Nt - 1

    t_string = T_string(nt);
    u_string = U_string(nt);

    % Find index of equal sign and extract time and velocity
    t = str2double(extractAfter(t_string, '='));
    u = str2double(extractAfter(u_string, '='));
    U_norm(nt, :) = [t * tc u / tc];

end

id = 1 : 10 : length(U_norm);
plot(U_norm(id, 1), U_norm(id, 2))
xlabel('$t$')
ylabel('$||\mathbf{u}||_2$')
xlim([0 300])
