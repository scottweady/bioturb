
cd scripts

for num = 1 : 8
    fprintf('Making figure %d...\n', num);
    run(sprintf('fig%d.m', num));
    print(gcf, sprintf('../fig%d', num), '-dpdf');
end

cd ..